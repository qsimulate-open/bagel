//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cis.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/response/cis.h>
#include <src/prop/multipole.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/davidson.h>

using namespace std;
using namespace bagel;

CIS::CIS(shared_ptr<const PTree> idata, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(idata, geom, ref), nstate_(idata->get<int>("nstate", 1)), nocc_(ref_->nocc()), nvirt_(ref_->nvirt()), maxiter_(idata->get<int>("maxiter", 20)),
    thresh_(idata->get<double>("thresh", 1.0e-6)), eig_(nocc_+nvirt_) {

  if (nocc_+nvirt_ != ref->coeff()->mdim())
    throw runtime_error("nocc + nvirt does not match the dimension of the coefficient");

  const MatView ocoeff = ref->coeff()->slice(0, nocc_);

  // compute half-transformed integrals
  half_ = geom_->df()->compute_half_transform(ocoeff);
  shared_ptr<const DFHalfDist> halfjj = half_->apply_JJ();

  shared_ptr<Matrix> fock = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, ocoeff, false/*dograd*/, true/*rhf*/);
  *fock = *ref->coeff() % *fock * *ref->coeff();
  fock->diagonalize(eig_);
  coeff_ = make_shared<Matrix>(*ref->coeff() * *fock);

  // re-compute half-transformed integrals
  half_ = geom_->df()->compute_half_transform(coeff_->slice(0, nocc_));
  fulljj_ = half_->compute_second_transform(coeff_->slice(0, nocc_))->apply_JJ();

}


void CIS::compute() {
  // initial guess
  {
    Matrix ebj(nvirt_, nocc_);
    for (int i = 0; i != nocc_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        ebj(j, i) = eig_[nocc_+j] - eig_[i];

    for (int i = 0; i != nstate_; ++i) {
      auto min = min_element(ebj.begin(), ebj.end());
      auto amp = make_shared<Matrix>(nvirt_, nocc_);
      *(amp->data() + distance(ebj.begin(), min)) += 1.0;
      amp_.push_back(amp);
      *min = 1.0e100;
    }
  }

  DavidsonDiag<Matrix> davidson(nstate_, maxiter_);

  const MatView ocoeff = coeff_->slice(0, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  vector<bool> conv(nstate_, false);

  cout << "  === CIS iteration ===" << endl << endl;
  Timer timer;

  for (int iter = 0; iter != maxiter_; ++iter) {
    vector<shared_ptr<const Matrix>> sigma;
    for (int ist = 0; ist != nstate_; ++ist) {
      if (!conv[ist]) {
        shared_ptr<Matrix> tmp = amp_[ist]->copy();
        // one body part
        for (int i = 0; i != nocc_; ++i)
          for (int j = 0; j != nvirt_; ++j)
            (*tmp)(j, i) *= eig_[nocc_+j] - eig_[i];
        const Matrix ovcoeff(vcoeff * *amp_[ist]);
        // J-type term
        Matrix one_occ(*geom_->df()->compute_Jop(half_, make_shared<Matrix>(*ovcoeff.transpose()*2.0), false) * ocoeff);
        // K-type term
        auto chalf = geom_->df()->compute_half_transform(ovcoeff);
        one_occ += *chalf->form_2index(fulljj_, -1.0);
        *tmp += vcoeff % one_occ;
        sigma.push_back(tmp);
      } else {
        sigma.push_back(nullptr);
      }
    }
    assert(amp_.size() == sigma.size());

    energy_ = davidson.compute(amp_, sigma);
    vector<shared_ptr<Matrix>> residual = davidson.residual();

    amp_.clear();
    for (int ist = 0; ist != nstate_; ++ist) {
      const double rms = residual[ist]->rms();
      cout << "      " << setw(3) << iter << (conv[ist] ? " * " : "   ")
           << setw(3) << ist << setw(15) << setprecision(8) << energy_[ist] << setw(15) << rms << setw(10) << setprecision(2) << timer.tick() << endl;
      if (rms > thresh_) {
        for (int i = 0; i != nocc_; ++i)
          for (int j = 0; j != nvirt_; ++j)
            residual[ist]->element(j, i) /= (eig_[nocc_+j] - eig_[i] - energy_[ist]);
        residual[ist]->scale(1.0/residual[ist]->norm());
        amp_.push_back(residual[ist]);
        conv[ist] = false;
      } else {
        amp_.push_back(nullptr);
        conv[ist] = true;
      }
    }
    if (nstate_ > 1)
      cout << endl;

    if (all_of(conv.begin(), conv.end(), [](const bool i) { return i; })) {
      cout << "    * CIS converged" << endl << endl;
      break;
    }
  }

  // setting CI amplitudes. Also compute transition dipole moments
  vector<shared_ptr<Matrix>> amp = davidson.civec();
  for (int ist = 0; ist != nstate_; ++ist) {
    amp_[ist] = amp[ist];
    auto transao = make_shared<Matrix>(vcoeff * *amp_[ist] ^ ocoeff);
    stringstream ss; ss << "Transition 0 -> " << ist;
    Dipole dipole(geom_, transao, ss.str());
    dipole.compute();
  }
}
