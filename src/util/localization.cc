//
// BAGEL - Parallel electron correlation program.
// Filename: localization.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/util/localization.h>
#include <src/scf/overlap.h>
#include <src/scf/coeff.h>
#include <src/util/jacobi.h>

using namespace bagel;
using namespace std;

/************************************************************************************
* Orthogonalize based on regions - follows the implementation in                    *
*   de Silva, Giebultowski, Korchowiec, PCCP 2011, 14, 546–552                      *
************************************************************************************/
RegionLocalization::RegionLocalization(shared_ptr<const Geometry> geom, shared_ptr<Matrix> density, vector<pair<int, int> > atom_bounds)
   : OrbitalLocalization(geom, density) {
  for (pair<int,int>& ibound : atom_bounds) {
    int start = geom_->offset(ibound.first).front();
    int end = (((ibound.second+1) >= (geom_->natom())) ? geom_->nbasis() : geom_->offset(ibound.second+1).front() );
    bounds_.push_back(make_pair(start, end));
  }
}

shared_ptr<Matrix> RegionLocalization::localize() {
  const int nbasis = geom_->nbasis();

  Overlap sqrt_S(geom_); sqrt_S.sqrt();

  // Symmetric orthogonalized density matrix
  shared_ptr<Matrix> ortho_density(new Matrix( sqrt_S * (*density_) * sqrt_S ));

  // transform will hold the eigenvectors of each block
  shared_ptr<Matrix> T(new Matrix(nbasis, nbasis));
  vector<double> eigenvalues;

  // Go through each region and diagonalize the blocks, the eigenvectors form T
  for(auto& iregion : bounds_) {
    const int start = iregion.first;
    const int fence = iregion.second;
    const int size_i = fence - start;

    shared_ptr<Matrix> D_i = ortho_density->get_submatrix(start, start, size_i, size_i);

    vector<double> eig_i(size_i);
    D_i->diagonalize(eig_i.data());

    eigenvalues.insert(eigenvalues.end(), eig_i.begin(), eig_i.end());
    T->copy_block(start, start, size_i, size_i, D_i->data());
  }

  *ortho_density = (*T) % (*ortho_density) * (*T);

  // Set up for Jacobi rotations. Maybe this is general enough to have its own utility?

  // U matrix will collect the transformations. It should start as the identity.
  shared_ptr<Matrix> U(new Matrix(nbasis, nbasis)); U->add_diag(1.0);

  // Classify each eigenvalue as occupied, mixed, or virtual. All of these classifications may need to be revisited at some point.
  vector<int> occupied, mixed, virt;
  {
    int imo = 0;
    for(auto eig_iter = eigenvalues.begin(); eig_iter != eigenvalues.end(); ++eig_iter, ++imo) {
      if (*eig_iter > 1.5) occupied.push_back(imo);
      else if (*eig_iter > 0.5) mixed.push_back(imo);
      else virt.push_back(imo);
    }
  }

  // Starting rotations. For now, let's see if one sweep suffices.
  {
    shared_ptr<JacobiDiag> jacobi(new JacobiDiag(ortho_density,U));

    for(int& iocc : occupied) {
      for(int& imixed : mixed) jacobi->rotate(iocc, imixed);
      for(int& ivirt : virt) jacobi->rotate(iocc, ivirt);
    }
  
    for(int& imixed : mixed) {
      for(int& ivirt : virt) jacobi->rotate(imixed, ivirt);
    }
  }

  // TODO: here would be the stage 2 jacobi sweeps. This would be necessary for covalently bound regions.

  Overlap S_inverse_half(geom_); S_inverse_half.inverse_half();

  // Reorder so that the occupied orbitals come first, separated by fragment
  shared_ptr<Matrix> tmp(new Matrix(S_inverse_half * (*T) * (*U)));

  shared_ptr<Matrix> out(new Matrix(nbasis, nbasis));

  {
    int imo = 0;
    for(int& iocc : occupied) {
      copy_n(tmp->element_ptr(0,iocc), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
    for(int& imixed : mixed) {
      copy_n(tmp->element_ptr(0,imixed), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
    for(int& ivirt : virt) {
      copy_n(tmp->element_ptr(0,ivirt), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
  }

  return out;
}

PMLocalization::PMLocalization(shared_ptr<const Geometry> geom, shared_ptr<const Coeff> coeff, const int norb) :
  OrbitalLocalization(geom, shared_ptr<Matrix>(new Matrix(*coeff))), coeff_(coeff), norb_(norb) {
  
  S_ = shared_ptr<Matrix>(new Overlap(geom));

  vector<vector<int> > offsets = geom->offsets();
  for(auto ioffset = offsets.begin(); ioffset != offsets.end(); ++ioffset) {
    const int start = ioffset->front();
    int end;
    if((ioffset+1) == offsets.end()) end = geom->nbasis();
    else end = (ioffset+1)->front();

    atom_bounds_.push_back(make_pair(start, end));
  }
}

shared_ptr<Matrix> PMLocalization::localize() {
  shared_ptr<JacobiPM> jacobi(new JacobiPM(density_, norb_, S_, atom_bounds_));

  cout << " === Starting Pipek-Mezey Localization ===" << endl << endl;
  cout << "Iteration          P              dP" << endl;

  double P = calc_P();

  cout << setw(3) << 0 << setw(16) << setprecision(10) << P << endl;

  for(int iter = 0; iter < 10; ++iter) {
    jacobi->sweep();

    double tmp_P = calc_P();
    cout << setw(3) << iter+1 << setw(16) << setprecision(10) << tmp_P
                              << setw(16) << setprecision(10) << tmp_P - P << endl;
    P = tmp_P;
  }

  return density_;
}

double PMLocalization::calc_P() const {
  double out = 0.0;

  const int nbasis = density_->ndim();
  for(int imo = 0; imo < norb_; ++imo) {
    for(auto& ibounds : atom_bounds_) {
      double QA = 0.0;
      for(int iao = ibounds.first; iao < ibounds.second; ++iao) {
        QA += density_->element(iao,imo) * ddot_(nbasis, density_->element_ptr(0,imo), 1, S_->element_ptr(0,iao), 1);
      }
      out += QA*QA;
    }
  }

  return out;
}
