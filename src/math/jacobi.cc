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


#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <src/math/matrix.h>
#include <src/math/jacobi.h>

using namespace std;
using namespace bagel;

void Jacobi_base::sweep() {
  for ( auto& isubsweep : sweeper_ ) {
    // These could all be done independently, in principle
    subsweep(isubsweep);
  }
}

void JacobiDiag::subsweep(vector<pair<int,int>> pairlist) {
  for (auto& ipair : pairlist) rotate(ipair.first, ipair.second);
}

void JacobiDiag::rotate(const int k, const int l) {
  // This is a hugely wasteful short term solution
  shared_ptr<Matrix> Alocal = A_->matrix();
  const double kl = Alocal->element(k,l);
  if (fabs(k) < numerical_zero__) return;

  const double kk = Alocal->element(k,k);
  const double ll = Alocal->element(l,l);

  const double beta = 0.5*(ll - kk)/kl;
  const double t = copysign(1.0,beta)/(fabs(beta) + sqrt(beta*beta + 1.0));
  const double c = 1.0/(sqrt(t*t + 1.0));
  const double s = c*t;
  const double rho = (1.0 - c)/s;

  Alocal->element(k,k) = kk - t * kl;
  Alocal->element(l,l) = ll + t * kl;

  Alocal->element(k,l) = 0.0;
  Alocal->element(l,k) = 0.0;

  // I'm afraid of overwriting data, thus copying some stuff
  unique_ptr<double[]> k_column(new double[nbasis_]);
  double* k_column_data = k_column.get();
  copy_n(Alocal->element_ptr(0,k), nbasis_, k_column_data);
  unique_ptr<double[]> l_column(new double[nbasis_]);
  double* l_column_data = l_column.get();
  copy_n(Alocal->element_ptr(0,l), nbasis_, l_column_data);

  for(int i = 0; i < nbasis_; ++i) {
    if (i == k || i == l) continue;

    const double ik = k_column_data[i];
    const double il = l_column_data[i];

    double new_ik = ik - s * (il + rho * ik);
    double new_il = il + s * (ik - rho * il);

    Alocal->element(i,k) = new_ik; Alocal->element(k,i) = new_ik;
    Alocal->element(i,l) = new_il; Alocal->element(l,i) = new_il;
  }

  A_ = Alocal->distmatrix();
  Q_->rotate(k, l, acos(c));
}

void JacobiPM::subsweep(vector<pair<int,int>> pairlist) {
  auto mos = make_shared<DistMatrix>(nbasis_, norb_);
  auto P_A = make_shared<DistMatrix>(norb_, norb_);

#ifdef HAVE_SCALAPACK
  pdgemm_("N", "N", nbasis_, norb_, nbasis_, 1.0, S_->local().get(), 1, 1, S_->desc().get(),
                                                 Q_->local().get(), 1, nstart_ + 1, Q_->desc().get(),
                                            0.0, mos->local().get(), 1, 1, mos->desc().get());
#else
  dgemm_("N", "N", nbasis_, norb_, nbasis_, 1.0, S_->data(), nbasis_, Q_->element_ptr(0, nstart_), nbasis_, 0.0, mos->data(), nbasis_);
#endif
  const int npairs = pairlist.size();

  vector<double> AA(npairs, 0.0);
  vector<double> BB(npairs, 0.0);

  for (auto& ibounds : atom_bounds_) {
    const int natombasis = ibounds.second - ibounds.first;

#ifdef HAVE_SCALAPACK
    pdgemm_("T", "N", norb_, norb_, natombasis, 1.0, mos->local().get(), ibounds.first + 1, 1, mos->desc().get(),
                                                         Q_->local().get(), ibounds.first + 1, nstart_ + 1, Q_->desc().get(),
                                                    0.0, P_A->local().get(), 1, 1, P_A->desc().get());

#else
    dgemm_("T", "N", norb_, norb_, natombasis, 1.0, mos->element_ptr(ibounds.first, 0), nbasis_,
                              Q_->element_ptr(ibounds.first, nstart_), nbasis_, 0.0, P_A->data(), norb_);
#endif

    shared_ptr<const Matrix> localPA = P_A->matrix(); // This should only be nocc x nocc in size

    for (int ipair = 0; ipair < npairs; ++ipair) {
      const int& kk = pairlist[ipair].first - nstart_;
      const int& ll = pairlist[ipair].second - nstart_;
      const double Qkl_A = 0.5 * (localPA->element( kk, ll ) + localPA->element( ll, kk ));
      const double Qkminusl_A = localPA->element(kk, kk) - localPA->element(ll, ll);

      AA[ipair] += Qkl_A*Qkl_A - 0.25*Qkminusl_A*Qkminusl_A;
      BB[ipair] += Qkl_A*Qkminusl_A;
    }
  }

  // This loop could be threaded
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
}
