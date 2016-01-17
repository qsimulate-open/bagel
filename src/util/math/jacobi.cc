//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: jacobi.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <functional>

#include <src/util/parallel/staticdist.h>
#include <src/util/math/matrix.h>
#include <src/util/math/jacobi.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;

void Jacobi_base::sweep() {
  for ( auto& isubsweep : sweeper_ ) {
    // These could all be done independently, in principle
    subsweep(isubsweep);
  }
}

void JacobiDiag::subsweep(vector<pair<int,int>>& pairlist) {
  for (auto& ipair : pairlist) rotate(ipair.first, ipair.second);
}

void JacobiDiag::rotate(const int k, const int l) {
  const double kl = A_->element(k,l);
  if (fabs(k) < numerical_zero__) return;

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

  Q_->rotate(k, l, acos(c));
}

void JacobiPM::subsweep(vector<pair<int,int>>& pairlist) {
  const int npairs = pairlist.size();

  // For parallelization across nodes
  StaticDist dist(npairs, mpi__->size());
  size_t pstart, pend;
  tie(pstart, pend) = dist.range(mpi__->rank());
  const size_t psize = pend - pstart;

  auto right = make_shared<Matrix>(nbasis_, 2*psize, true);
  for (size_t ip = 0; ip < psize; ++ip) {
    copy_n(SQ_->element_ptr(0, pairlist[ip + pstart].first), nbasis_, right->element_ptr(0, 2*ip));
    copy_n(SQ_->element_ptr(0, pairlist[ip + pstart].second), nbasis_, right->element_ptr(0, 2*ip+1));
  }

  shared_ptr<Matrix> left = ( lowdin_ ? right : make_shared<Matrix>(nbasis_, 2*psize, true) );
  if (!lowdin_) {
    for (size_t ip = 0; ip < psize; ++ip) {
      copy_n(Q_->element_ptr(0, pairlist[ip + pstart].first), nbasis_, left->element_ptr(0, 2*ip));
      copy_n(Q_->element_ptr(0, pairlist[ip + pstart].second), nbasis_, left->element_ptr(0, 2*ip+1));
    }
  }

  Matrix P_A(2*psize, 2*psize, true);

  vector<double> AA(npairs, 0.0);
  vector<double> BB(npairs, 0.0);

  for (auto& ibounds : atom_bounds_) {
    const int natombasis = ibounds.second - ibounds.first;
    const int boundstart = ibounds.first;

    for (int ip = 0; ip < psize; ++ip) {
      dgemm_("T", "N", 2, 2, natombasis, 1.0, left->element_ptr(boundstart, 2*ip), left->ndim(),
          right->element_ptr(boundstart, 2*ip), right->ndim(), 0.0, P_A.element_ptr(2*ip,2*ip), P_A.ndim());

      const int kk = 2*ip;
      const int ll = 2*ip+1;

      const double Qkl_A = 0.5 * (P_A(kk,ll) + P_A(ll,kk));
      const double Qkminusl_A = P_A(kk,kk) - P_A(ll,ll);

      AA[ip + pstart] += Qkl_A*Qkl_A - 0.25*Qkminusl_A*Qkminusl_A;
      BB[ip + pstart] += Qkl_A*Qkminusl_A;
    }
  }

  mpi__->allreduce(AA.data(), AA.size());
  mpi__->allreduce(BB.data(), BB.size());

  vector<tuple<int, int, double>> rotations;
  for (int ipair = 0; ipair < npairs; ++ipair) {
    const double Akl = AA[ipair];
    const double Bkl = BB[ipair];

    const int kk = pairlist[ipair].first;
    const int ll = pairlist[ipair].second;

    if( fabs(Bkl) < numerical_zero__ && Akl > 0.0 ) continue;

    double gamma = copysign(0.25, Bkl) * acos( -Akl/hypot(Akl,Bkl) );

    rotations.emplace_back(kk, ll, gamma);
  }

  Q_->rotate(rotations);
  SQ_->rotate(rotations);
}
