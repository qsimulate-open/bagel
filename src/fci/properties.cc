//
// BAGEL - Parallel electron correlation program.
// Filename: properties.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/fci/properties.h>
#include <src/molecule/dipolematrix.h>
#include <src/fci/prop1etask.h>

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
  for(int i = 0; i < 3; ++i) compressed_dipoles_[i] = unique_ptr<double[]>(new double[sizeij_]);

  for (int imu = 0; imu < 3; ++imu) {
    for (int i = 0, ij = 0; i < norb_; ++i) {
      for (int j = 0; j < norb_; ++j, ++ij) {
        compressed_dipoles_[imu][ij] = dipole_mo_[imu]->element(i + nocc_, j + nocc_);
      }
    }
  }
}

void CIDipole::compute(std::shared_ptr<const Dvec> ccvec) {
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  auto sigma = make_shared<Dvec>(det, nstates);

  shared_ptr<const Determinants> det_trans = det->transpose();
  shared_ptr<const Dvec> cc_trans = ccvec->spinflip(det_trans);
  auto sg_trans = make_shared<Dvec>(det, nstates);

  const int sizeij = sizeij_;
  const int la = ccvec->lena();
  const int lb = ccvec->lenb();

  vector<Prop1eTask> tasks;
  tasks.reserve(nstates * (la + lb));

  for (int imu = 0; imu < 3; ++imu) {
    sigma->zero();
    sg_trans->zero();

    for (int istate = 0; istate < nstates; ++istate) {
      double *target = sigma->data(istate)->data();
      for (auto aiter = det->stringa().begin(); aiter != det->stringa().end(); ++aiter, target+=lb)
        tasks.emplace_back(ccvec->data(istate), *aiter, target, compressed_dipoles_[imu].get());

      target = sg_trans->data(istate)->data();
      for (auto aiter = det_trans->stringa().begin(); aiter != det_trans->stringa().end(); ++aiter, target+=la)
        tasks.emplace_back(ccvec->data(istate), *aiter, target, compressed_dipoles_[imu].get());
    }

    TaskQueue<Prop1eTask> tq(tasks);
    tq.compute(resources__->max_num_threads());

    auto tmp = make_shared<Matrix>(nstates, nstates);
    for (int j = 0; j < nstates; ++j) {
      for (int i = 0; i < nstates; ++i) {
        // This saves an extra transpose at the expense of more ddot... which is better?
        tmp->element(i,j) = sigma->data(i)->ddot(*ccvec->data(j)) + sg_trans->data(i)->ddot(*cc_trans->data(j));
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
