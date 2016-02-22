//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: properties.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/ci/fci/properties.h>
#include <src/ci/fci/prop1etask.h>
#include <src/mat1e/dipolematrix.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIProperties)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Prop1e)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIDipole)

using namespace std;
using namespace bagel;

// CIDipole start
void CIDipole::init(const int nstart, const int nfence) {
  DipoleMatrix dipole(geom_);

  for(int i = 0; i < dipole.Nblocks(); ++i) {
    dipole_mo_[i] = make_shared<Matrix>((*coeff_) % (*dipole.data(i)) * (*coeff_));
  }

  core_dipole_ = {{0.0, 0.0, 0.0}} ;
  for(int i = 0; i < dipole.Nblocks(); ++i) {
    for(int j = 0; j < nocc_; ++j) {
      core_dipole_[i] += 2.0 * dipole_mo_[i]->element(j,j);
    }
  }

  sizeij_ = norb_ * norb_;

  for (int imu = 0; imu < 3; ++imu)
    compressed_dipoles_[imu] = dipole_mo_[imu]->get_submatrix(nocc_, nocc_, norb_, norb_);
}

void CIDipole::compute(shared_ptr<const Dvec> ccvec) {
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  auto sigma = make_shared<Dvec>(det, nstates);

  shared_ptr<const Determinants> det_trans = det->transpose();
  shared_ptr<const CASDvec> cc = make_shared<CASDvec>(ccvec->dvec());
  shared_ptr<const CASDvec> cc_trans = cc->transpose(det_trans);
  auto sg_trans = make_shared<Dvec>(det, nstates);

  const int la = ccvec->lena();
  const int lb = ccvec->lenb();

  TaskQueue<Prop1eTask> tasks(nstates * (la + lb));

  for (int imu = 0; imu < 3; ++imu) {
    sigma->zero();
    sg_trans->zero();

    for (int istate = 0; istate < nstates; ++istate) {
      double *target = sigma->data(istate)->data();
      for (auto& a : det->string_bits_a()) {
        tasks.emplace_back(ccvec->data(istate), a, target, compressed_dipoles_[imu]->data());
        target += lb;
      }

      target = sg_trans->data(istate)->data();
      for (auto& a : det_trans->string_bits_a()) {
        tasks.emplace_back(ccvec->data(istate), a, target, compressed_dipoles_[imu]->data());
        target += la;
      }
    }

    tasks.compute();

    auto tmp = make_shared<Matrix>(nstates, nstates);
    for (int j = 0; j < nstates; ++j) {
      for (int i = 0; i < nstates; ++i) {
        // This saves an extra transpose at the expense of more ddot... which is better?
        tmp->element(i,j) = sigma->data(i)->dot_product(*ccvec->data(j)) + sg_trans->data(i)->dot_product(*cc_trans->data(j));
      }
    }
    tmp->add_diag(core_dipole_[imu]);
    // Maybe scale into different units?
    dipole_matrices_[imu] = tmp;
  }
}

void CIDipole::print() const {
  const int nstates = dipole_matrices_[0]->ndim();

  const string indent("   ");

  cout << endl << indent << " --- Dipole Moments ---" << endl;
  cout << indent << " mu_x    |0>";
  for (int istate = 1; istate < nstates; ++istate) cout << "         |" << istate << ">";
  cout << endl;
  for (int istate = 0; istate < nstates; ++istate) {
    cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nstates; ++jstate) {
      cout << setw(12) << setprecision(6) << dipole_matrices_[0]->element(jstate, istate);
    }
    cout << endl;
  }
  cout << endl;

  cout << indent << " mu_y    |0>";
  for (int istate = 1; istate < nstates; ++istate) cout << "         |" << istate << ">";
  cout << endl;
  for (int istate = 0; istate < nstates; ++istate) {
    cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nstates; ++jstate) {
      cout << setw(12) << setprecision(6) << dipole_matrices_[1]->element(jstate, istate);
    }
    cout << endl;
  }
  cout << endl;

  cout << indent << " mu_z    |0>";
  for (int istate = 1; istate < nstates; ++istate) cout << "         |" << istate << ">";
  cout << endl;
  for (int istate = 0; istate < nstates; ++istate) {
    cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nstates; ++jstate) {
      cout << setw(12) << setprecision(6) << dipole_matrices_[2]->element(jstate, istate);
    }
    cout << endl;
  }
  cout << endl;
}
