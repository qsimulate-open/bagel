//
// BAGEL - Parallel electron correlation program.
// Filename: dirac_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/util/constants.h>
#include <src/london/dirac_london.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <src/math/diis.h>
#include <src/rel/dfock.h>
#include <src/rel/relreference.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Dirac_London)

Dirac_London::Dirac_London(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom,
             const shared_ptr<const Reference> re) : Method(idata, geom, re) {
  gaunt_ = idata->get<bool>("gaunt", false);
  breit_ = idata->get<bool>("breit", gaunt_);
  robust_ = idata->get<bool>("robust", false);

  // when computing gradient, we store half-transform integrals
  do_grad_ = idata->get<bool>("gradient", false);

  geom_ = geom->relativistic(gaunt_);
  common_init(idata);
}


void Dirac_London::common_init(const shared_ptr<const PTree> idata) {
  Timer init_time;
  cout << "  *** Dirac HF ***" << endl << endl;

  // reading input keywords
  max_iter_ = idata->get<int>("maxiter", 100);
  max_iter_ = idata->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata->get<int>("diis_start", 1);
  thresh_scf_ = idata->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata->get<double>("thresh_scf", thresh_scf_);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  ncharge_ = idata->get<int>("charge", 0);
  nele_ = geom_->nele()-ncharge_;

  hcore_ = make_shared<const RelHcore_London>(geom_);
  init_time.tick_print("1-electron integrals");

  overlap_ = make_shared<const RelOverlap_London>(geom_);
  s12_ = overlap_->tildex(thresh_overlap_);

  nneg_ = s12_->mdim()/2;
  assert(s12_->mdim() % 2 == 0);

  if (breit_ && !gaunt_) throw runtime_error("Breit cannot be turned on if Gaunt is off");
}


void Dirac_London::compute() {
  Timer scftime;
  string indent = "  ";

  shared_ptr<const DistZMatrix> hcore = hcore_->distmatrix();
  shared_ptr<const DistZMatrix> distovl = overlap_->distmatrix();
  shared_ptr<const DistZMatrix> s12 = s12_->distmatrix();
  eig_ = VectorB(hcore->ndim());

  // making initial guess
  shared_ptr<const DistZMatrix> coeff = initial_guess(s12, hcore);
  shared_ptr<const DistZMatrix> aodensity = coeff->form_density_rhf(nele_, nneg_);

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl << endl;
  cout << indent << "    * DIIS with orbital gradients will be used." << endl << endl;
  scftime.tick_print("SCF startup");
  cout << endl;
  cout << indent << "=== Dirac RHF iteration (" + geom_->basisfile() + ", RMB) ===" << endl << indent << endl;

  DIIS<DistZMatrix, ZMatrix> diis(5);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer ptime(1);

    auto fock = make_shared<DFock>(geom_, hcore_, coeff->matrix()->slice_copy(nneg_, nele_+nneg_), gaunt_, breit_, do_grad_, robust_);

// TODO I have a feeling that the code should not need this, but sometimes there are slight errors. still looking on it.
#if 0
    fock->hermite();
#endif
    // distribute
    shared_ptr<const DistZMatrix> distfock = fock->distmatrix();

    // compute energy here
    const complex<double> prod = aodensity->dot_product(*hcore+*distfock); // identical to Tr(D^+ F)
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << prod.imag();
//    throw runtime_error(ss.str());
      cout << ss.str() << endl;
    }
    energy_ = 0.5*prod.real() + geom_->nuclear_repulsion();

    auto error_vector = make_shared<const DistZMatrix>(*distfock**aodensity**distovl - *distovl**aodensity**distfock);
    const double error = error_vector->rms();

    ptime.tick_print("Fock build");
    cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_
         << "   " << setw(17) << error << setw(15) << setprecision(2) << scftime.tick() << endl;

    if (error < thresh_scf_ && iter > 0) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      // when computing gradient, we store half-transform integrals to avoid recomputation
      if (do_grad_) {
        throw std::logic_error("Gradient integrals have not been implemented with London orbitals yet.");
        half_ = fock->half();
      }
      break;
    } else if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in SCF." << endl << endl;
      throw runtime_error("Max iteration reached in Dirac--Fock SCF");
    }

    if (iter >= diis_start_) {
      distfock = diis.extrapolate({distfock, error_vector});
      ptime.tick_print("DIIS");
    }

    DistZMatrix intermediate(*coeff % *distfock * *coeff);
    intermediate.diagonalize(eig_);
    coeff = make_shared<DistZMatrix>(*coeff * intermediate);

    aodensity = coeff->form_density_rhf(nele_, nneg_);

  }

  coeff_ = coeff->matrix();
}


//Print non dirac sea eigenvalues
void Dirac_London::print_eig() const {
  const int n = geom_->nbasis();
  for (int i = 0*n; i != 4*n; ++i) cout << setprecision(10) << setw(15) << eig_[i] <<  endl;
}


shared_ptr<const Reference> Dirac_London::conv_to_ref() const {
  // we store only positive state coefficients
  const size_t npos = coeff_->mdim() - nneg_;
  auto out = make_shared<RelReference>(geom_, coeff_, energy_, nneg_, nele_, npos-nele_, gaunt_, breit_, true);
  vector<double> eigp(eig_.begin()+nneg_, eig_.end());
  vector<double> eigm(eig_.begin(), eig_.begin()+nneg_);
  VectorB eig(eig_.size());
  copy(eigp.begin(), eigp.end(), eig.begin());
  copy(eigm.begin(), eigm.end(), eig.begin()+eigp.size());
  out->set_eig(eig);
  return out;
}


shared_ptr<const DistZMatrix> Dirac_London::initial_guess(const shared_ptr<const DistZMatrix> s12, const shared_ptr<const DistZMatrix> hcore) const {
  const int n = geom_->nbasis();
  VectorB eig(hcore->ndim());

  shared_ptr<const DistZMatrix> coeff;
  if (!ref_) {
    // No reference; starting from hcore
    DistZMatrix interm = *s12 % *hcore * *s12;
    interm.diagonalize(eig);
    coeff = make_shared<const DistZMatrix>(*s12 * interm);

  } else if (dynamic_pointer_cast<const RelReference>(ref_)) {
    auto relref = dynamic_pointer_cast<const RelReference>(ref_);

    if (relref->rel()) {
      // Relativistic, GIAO-based reference
      shared_ptr<ZMatrix> fock = make_shared<DFock>(geom_, hcore_, relref->relcoeff()->slice_copy(0, nele_), gaunt_, breit_, /*store_half*/false, robust_);
      DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
      interm.diagonalize(eig);
      coeff = make_shared<const DistZMatrix>(*s12 * interm);

    } else if (!relref->rel()) {
      // Non-relativistic, GIAO-based reference
      const int nocc = ref_->nocc();
      shared_ptr<ZMatrix> fock;
      assert(nocc*2 == nele_);
      auto ocoeff = make_shared<ZMatrix>(n*4, 2*nocc);
      ocoeff->add_block(1.0, 0,    0, n, nocc, relref->relcoeff()->slice(0,nocc));
      ocoeff->add_block(1.0, n, nocc, n, nocc, relref->relcoeff()->slice(0,nocc));
      fock = make_shared<DFock>(geom_, hcore_, ocoeff, gaunt_, breit_, /*store_half*/false, robust_);
      DistZMatrix interm = *s12 % *fock->distmatrix() * *s12;
      interm.diagonalize(eig);
      coeff = make_shared<const DistZMatrix>(*s12 * interm);
    }
  } else {
    throw logic_error("Invalid Reference passed to Dirac_London");
  }
  return coeff;
}
