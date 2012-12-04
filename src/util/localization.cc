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
#include <src/util/f77.h>
#include <src/util/localization.h>
#include <src/scf/overlap.h>
#include <src/util/constants.h>

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

    Matrix D_i(size_i, size_i);
    D_i.copy_block(0,0,size_i,size_i, ortho_density->get_block(start, start, size_i, size_i));

    vector<double> eig_i(size_i);
    D_i.diagonalize(eig_i.data());

    eigenvalues.insert(eigenvalues.end(), eig_i.begin(), eig_i.end());
    T->copy_block(start, start, size_i, size_i, D_i.data());
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
    for(int& iocc : occupied) {
      for(int& imixed : mixed) jacobi(iocc, imixed, ortho_density, U);
      for(int& ivirt : virt) jacobi(iocc, ivirt, ortho_density, U);
    }
  
    for(int& imixed : mixed) {
      for(int& ivirt : virt) jacobi(imixed, ivirt, ortho_density, U);
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

  cout << "The diagonal components of the transformed DM are: " << endl;
  for (int i = 0; i < nbasis; ++i) {
    cout << ortho_density->element(i,i);
    for (int j = 0; j < 7 && i < nbasis; ++j, ++i) {
      cout << ",   " << ortho_density->element(i,i);
    }
    cout << endl;
  }

  return out;
}

// For better or for worse, implementation following the description on wikipedia :)
void OrbitalLocalization::jacobi(const int k, const int l, shared_ptr<Matrix> A, shared_ptr<Matrix> Q) {
  const double kl = A->element(k,l);
  if (k < numerical_zero__) return;
  const double kk = A->element(k,k);
  const double ll = A->element(l,l);

  const bool vectors = (Q != nullptr);

  const double beta = 0.5*(ll - kk)/kl;
  const double t = copysign(1.0,beta)/(abs(beta) + sqrt(beta*beta + 1.0));
  const double c = 1.0/(sqrt(t*t + 1.0));
  const double s = c*t;
  const double rho = (1.0 - c)/s;
  
  A->element(k,k) = kk - t * kl;
  A->element(l,l) = ll + t * kl;

  A->element(k,l) = 0.0;
  A->element(l,k) = 0.0;

  const int nbasis = A->ndim();

  // I'm afraid of overwriting data, thus copying some stuff
  unique_ptr<double[]> k_column(new double[nbasis]);
  double* k_column_data = k_column.get();
  copy_n(A->element_ptr(0,k), nbasis, k_column_data);
  unique_ptr<double[]> l_column(new double[nbasis]);
  double* l_column_data = l_column.get();
  copy_n(A->element_ptr(0,l), nbasis, l_column_data);
  
  for(int i = 0; i < nbasis; ++i) {
    if (i == k || i == l) continue;

    const double ik = k_column_data[i];
    const double il = l_column_data[i];

    double new_ik = ik - s * (il + rho * ik);
    double new_il = il + s * (ik - rho * il);

    A->element(i,k) = new_ik; A->element(k,i) = new_ik;
    A->element(i,l) = new_il; A->element(l,i) = new_il;
  }

  if(vectors) drot_(nbasis, Q->element_ptr(0,k), 1, Q->element_ptr(0,l), 1, c, -s);
}
