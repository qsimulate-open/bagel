//
// BAGEL - Parallel electron correlation program.
// Filename: jacobi.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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
#include <src/wfn/geometry.h>
#include <src/util/matrix.h>
#include <src/util/jacobi.h>

using namespace std;
using namespace bagel;

void Jacobi_base::sweep() {
  for(int i = nstart_; i < (norb_ + nstart_); ++i) {
    for(int j = nstart_; j < i; ++j) {
      rotate(i,j);
    }
  }
}

void JacobiDiag::rotate(const int k, const int l) {
  const double kl = A_->element(k,l);
  if (abs(k) < numerical_zero__) return;

  const double kk = A_->element(k,k);
  const double ll = A_->element(l,l);

  const double beta = 0.5*(ll - kk)/kl;
  const double t = copysign(1.0,beta)/(fabs(beta) + sqrt(beta*beta + 1.0));
  const double c = 1.0/(sqrt(t*t + 1.0));
  const double s = c*t;
  const double rho = (1.0 - c)/s;

  A_->element(k,k) = kk - t * kl;
  A_->element(l,l) = ll + t * kl;

  A_->element(k,l) = 0.0;
  A_->element(l,k) = 0.0;

  // I'm afraid of overwriting data, thus copying some stuff
  unique_ptr<double[]> k_column(new double[nbasis_]);
  double* k_column_data = k_column.get();
  copy_n(A_->element_ptr(0,k), nbasis_, k_column_data);
  unique_ptr<double[]> l_column(new double[nbasis_]);
  double* l_column_data = l_column.get();
  copy_n(A_->element_ptr(0,l), nbasis_, l_column_data);

  for(int i = 0; i < nbasis_; ++i) {
    if (i == k || i == l) continue;

    const double ik = k_column_data[i];
    const double il = l_column_data[i];

    double new_ik = ik - s * (il + rho * ik);
    double new_il = il + s * (ik - rho * il);

    A_->element(i,k) = new_ik; A_->element(k,i) = new_ik;
    A_->element(i,l) = new_il; A_->element(l,i) = new_il;
  }

  drot_(nbasis_, Q_->element_ptr(0,k), 1, Q_->element_ptr(0,l), 1, c, -s);
}

void JacobiPM::rotate(const int k, const int l) {
  // There may be a more efficient way to do this, but this works for now

  vector<double> Qkl_A, Qkminusl_A;
  for(auto& ibounds : atom_bounds_) {
    double value_kl = 0.0, value_k = 0.0, value_l = 0.0;

    for(int i = ibounds.first; i < ibounds.second; ++i) {
      value_kl += Q_->element(i,l) * ddot_(nbasis_, Q_->element_ptr(0,k), 1, S_->element_ptr(0,i), 1);
      value_kl += Q_->element(i,k) * ddot_(nbasis_, Q_->element_ptr(0,l), 1, S_->element_ptr(0,i), 1);

      value_k += Q_->element(i,k) * ddot_(nbasis_, Q_->element_ptr(0,k), 1, S_->element_ptr(0,i), 1);
      value_l += Q_->element(i,l) * ddot_(nbasis_, Q_->element_ptr(0,l), 1, S_->element_ptr(0,i), 1);
    }

    Qkl_A.push_back(0.5*value_kl);
    Qkminusl_A.push_back(value_k-value_l);
  }

  const int natoms = atom_bounds_.size();

  double Ast = ddot_(natoms, Qkl_A.data(), 1, Qkl_A.data(), 1) - 0.25 * ddot_(natoms, Qkminusl_A.data(), 1, Qkminusl_A.data(), 1);
  double Bst = ddot_(natoms, Qkl_A.data(), 1, Qkminusl_A.data(), 1);

  if( fabs(Bst) < numerical_zero__ && Ast > 0.0 ) return;

  double gamma = copysign(0.25, Bst) * acos( -Ast/hypot(Ast,Bst) );

  drot_(nbasis_, Q_->element_ptr(0,k), 1, Q_->element_ptr(0,l), 1, cos(gamma), sin(gamma));
}
