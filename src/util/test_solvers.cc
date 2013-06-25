//
// BAGEL - Parallel electron correlation program.
// Filename: test_solvers.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include <src/scf/matrix1e.h>
#include <src/util/linear.h>
#include <src/util/linearRM.h>
#include <src/util/aughess.h>
#include <src/util/davidson.h>
#include <src/util/zdavidson.h>
#include <src/util/bfgs.h>
#include <src/scf/scf.h>
#include <iostream>

using namespace std;
using namespace bagel;

void test_solvers(shared_ptr<Geometry> geom_) {
{
 cout << " Testing solvers." << endl;
 shared_ptr<Matrix> target(new Matrix(geom_->nbasis(), geom_->nbasis()));

 assert(target->ndim() == target->mdim());
 const size_t n = target->ndim();

 // source term
 for (int i = 0; i != n*n; ++i) {
   target->data(i) = static_cast<double>(rand()) / RAND_MAX;
 }
 // matrix.
 unique_ptr<double[]> hess(new double[n*n*n*n]);
 shared_ptr<Matrix> diag(new Matrix(geom_->nbasis(), geom_->nbasis()));
 for (int i = 0; i != n*n*n*n; ++i) {
   hess[i] = static_cast<double>(rand()) / RAND_MAX;
 }
 for (int i = 0; i != n*n; ++i) {
   for (int j = 0; j != n*n; ++j) {
     hess[j+n*n*i] =  hess[i+n*n*j] = 0.5*(hess[j+n*n*i] + hess[i+n*n*j]);
   }
   hess[i+n*n*i] += 2.0*i;
   diag->data(i) = hess[i+n*n*i];
 }

 const double tiny = 1.0e-20;

#if 1
  // testing Davidson -- checked.
  cout << "  testing Davidson class" << endl;
  DavidsonDiag<Matrix> davidson(1,n*n);
  shared_ptr<Matrix> prev(new Matrix(n, n));
  prev->element(0,0) = 1.0;

  for (int i = 0; i != n*n; ++i) {
    shared_ptr<Matrix> start(new Matrix(*prev));
    davidson.orthog(start);
    shared_ptr<Matrix> res = start->clone();
    dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);

    shared_ptr<const Matrix> ss(new Matrix(*start));
    shared_ptr<const Matrix> rr(new Matrix(*res));
    const double energy = davidson.compute(ss, rr);
    shared_ptr<Matrix> residual = davidson.residual().front();

    cout << "davidson " << setw(20) << setprecision(10) << fixed << ::pow(residual->norm(),2.0) << " " << setw(20) << energy << endl;
    if (::pow(residual->norm(),2.0) < tiny) break;

    for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
    prev = residual;
  }
}


{
  cout << "Testing Mike's solver." << endl;
  shared_ptr<Matrix> target(new Matrix(geom_->nbasis(), geom_->nbasis()));
  assert(target->ndim() == target->mdim());
  const size_t n = target->ndim();
  // source term
  for (int i = 0; i != n*n; ++i) {
    target->data(i) = static_cast<double>(rand()) / RAND_MAX;
  }
  // matrix.
  unique_ptr<complex<double>[]> hess(new complex<double>[n*n*n*n]);
  shared_ptr<ZMatrix> diag(new ZMatrix(geom_->nbasis(), geom_->nbasis()));
  for (int i = 0; i != n*n*n*n; ++i) {
    double a = static_cast<double>(rand()) / RAND_MAX;
    double b = static_cast<double>(rand()) / RAND_MAX;
    hess[i] = complex<double>(a, b);
  }
  for (int i = 0; i != n*n; ++i) {
    for (int j = 0; j != n*n; ++j) {
      if (i==j) {
        hess[j+n*n*i] = complex<double>(hess[i+n*n*j].real(), 0);
        //std::cout << "diagonal" << i << ", " << j << setw(5) << hess[j+n*n*i] << ", " << hess[i+n*n*j] << std::endl;
      }
      else {
        hess[j+n*n*i] = complex<double>(hess[i+n*n*j].real(), -hess[i+n*n*j].imag());
        //std::cout << "off-diagonal" << i << ", " << j << setw(5) << hess[j+n*n*i] << ", " << hess[i+n*n*j] << std::endl;
      }
    }
    hess[i+n*n*i] += 2.0*i;
    diag->data(i) = hess[i+n*n*i];
  }

  const double tiny = 1.0e-20;

  // testing Davidson -- checked.

  cout << "  testing ZDavidson class" << endl;
  ZDavidsonDiag<ZMatrix> zdavidson(1,n*n);
  shared_ptr<ZMatrix> prev(new ZMatrix(n, n));
  prev->element(0,0) = complex<double>(1.0, -1.0);

  for (int i = 0; i != n*n; ++i) {
    shared_ptr<ZMatrix> start(new ZMatrix(*prev));
    zdavidson.orthog(start);
    shared_ptr<ZMatrix> res = start->clone();
    zgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);

    shared_ptr<const ZMatrix> ss(new ZMatrix(*start));
    shared_ptr<const ZMatrix> rr(new ZMatrix(*res));
    const double energy = zdavidson.compute(ss, rr);
    shared_ptr<ZMatrix> residual = zdavidson.residual().front();

    cout << "davidson " << setw(20) << setprecision(10) << fixed << ::pow(residual->norm(),2.0) << " " << setw(20) << energy << endl;
    if (::pow(residual->norm(),2.0) < tiny) break;

    for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
    prev = residual;
  }
}

#endif
}
#if 0
  // testing Linear
  {
    cout << "  testing Linear class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));
    Linear<Matrix> linear(n*n, tmp);

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(geom_->nbasis(), geom_->nbasis()));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      linear.orthog(prev);
      shared_ptr<Matrix> res = prev->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, prev->data(), 1, 0.0, res->data(), 1);

      shared_ptr<Matrix> residual = linear.compute_residual(prev, res);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      for (int i = 0; i != diag->size(); ++i) residual->data(i) /= diag->data(i);
      prev = residual;
    }
    linear.civec()->print();
  }

  // testing Linear2
  {
    cout << "  testing Linear2 class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));
    LinearRM<Matrix> linear(n*n, tmp);

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(geom_->nbasis(), geom_->nbasis()));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      linear.orthog(prev);
      shared_ptr<Matrix> res = prev->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, prev->data(), 1, 0.0, res->data(), 1);

      shared_ptr<Matrix> residual = linear.compute_residual(prev, res);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      for (int i = 0; i != diag->size(); ++i) residual->data(i) /= diag->data(i);
      prev = residual;
    }
    linear.civec()->print();
  }

#if 0
  // TODO to be checked!!!!!!!!!!!!
  // testing AugHess
  {
    cout << "  testing AugHess class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));
    AugHess<Matrix> linear(n*n+1, tmp);

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(geom_->nbasis(), geom_->nbasis()));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix> start(new Matrix(*prev));
      linear.orthog(start);

      shared_ptr<Matrix> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);

      shared_ptr<Matrix> residual = linear.compute_residual(start, res);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      cout << "lamba " << linear.eig() << endl;

      for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
      prev = residual;
    }
    linear.civec()->print();
  }
#endif

#if 0
  // checked.
  // testing straight-line quasi-newton.
  {
    cout << "  testing Straight quasi-Newton class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(n, n));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix> start(new Matrix(*prev));
      shared_ptr<Matrix> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);
      shared_ptr<Matrix> residual(new Matrix(*res - *tmp));
      for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      *prev -= *residual;
      if (::pow(residual->norm(),2.0) < tiny) break;
    }
    prev->print();
  }
#endif

#if 0
  // testing BFGS update
  {
    cout << "  testing BFGS class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(n, n));
    prev->element(0,0) = 1.0;

    BFGS<Matrix> bfgs(diag);

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix> start(new Matrix(*prev));
      shared_ptr<Matrix> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);
      shared_ptr<Matrix> residual0(new Matrix(*res - *tmp));

      shared_ptr<Matrix> residual = bfgs.extrapolate(residual0, prev);

      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      *prev -= *residual;
      if (::pow(residual->norm(),2.0) < tiny) break;
    }
    prev->print();
  }

  // testing Linear with BFGS updates
  {
    cout << "  testing Linear class" << endl;
    shared_ptr<Matrix> tmp(new Matrix(*target));
    Linear<Matrix> linear(n*n, tmp);

    // start with target/denom
    shared_ptr<Matrix> prev(new Matrix(n, n));
    prev->element(0,0) = 1.0;

    BFGS<Matrix> bfgs(diag);

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix> start(new Matrix(*prev));
      linear.orthog(start);
      shared_ptr<Matrix> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1);

      shared_ptr<Matrix> residual0 = linear.compute_residual(start, res);

      shared_ptr<Matrix> residual = bfgs.extrapolate(residual0, linear.civec());
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      prev = residual;
    }
    linear.civec()->print();
  }
#endif

#if 0
  // reference
  shared_ptr<Matrix> answer(new Matrix(*target));
  unique_ptr<int[]> ipiv(new int[n*n]);
  int info;
  dgesv_(n*n, 1, hess.get(), n*n, ipiv.get(), answer->data(), n*n, info);
  answer->print();
#endif
}

#else
static int a = 0; // just to have something in this file
#endif
