//
// Author : Toru Shiozaki
// Date   : April 2012
//

#include <src/scf/matrix1e.h>
#include <src/util/linear.h>
#include <src/util/aughess.h>
#include <src/util/davidson.h>
#include <src/util/bfgs.h>
#include <src/scf/scf.h>
#include <iostream>

using namespace std;

void test_solvers(shared_ptr<Geometry> geom_) {
  cout << " Testing solvers." << endl;
  shared_ptr<Matrix1e> target(new Matrix1e(geom_));

  assert(target->ndim() == target->mdim());
  const size_t n = target->ndim();

  // source term
  for (int i = 0; i != n*n; ++i) {
    target->data(i) = static_cast<double>(rand()) / RAND_MAX;
  }
  // matrix.
  unique_ptr<double[]> hess(new double[n*n*n*n]);
  shared_ptr<Matrix1e> diag(new Matrix1e(geom_));
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
  {
    cout << "  testing Davidson class" << endl;
    DavidsonDiag<Matrix1e> davidson(1,n*n);
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n*n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      davidson.orthog(start);
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 

      shared_ptr<const Matrix1e> ss(new Matrix1e(*start)); 
      shared_ptr<const Matrix1e> rr(new Matrix1e(*res)); 
      const double energy = davidson.compute(ss, rr);
      shared_ptr<Matrix1e> residual = davidson.residual().front();

      cout << "davidson " << setw(20) << setprecision(10) << fixed << ::pow(residual->norm(),2.0) << " " << setw(20) << energy << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
      prev = residual;
    }
  }
#endif

  // testing Linear
  {
    cout << "  testing Linear class" << endl;
    shared_ptr<Matrix1e> tmp(new Matrix1e(*target));
    Linear<Matrix1e> linear(n*n, tmp);

    // start with target/denom
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      linear.orthog(start);
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 

      shared_ptr<Matrix1e> residual = linear.compute_residual(start, res);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
      prev = residual;
    }
    linear.civec()->print();
  }

#if 0
  // TODO to be checked!!!!!!!!!!!!
  // testing AugHess
  {
    cout << "  testing AugHess class" << endl;
    shared_ptr<Matrix1e> tmp(new Matrix1e(*target));
    AugHess<Matrix1e> linear(n*n+1, tmp);

    // start with target/denom
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      linear.orthog(start);
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 

      shared_ptr<Matrix1e> residual = linear.compute_residual(start, res);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

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
    shared_ptr<Matrix1e> tmp(new Matrix1e(*target));

    // start with target/denom
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 
      shared_ptr<Matrix1e> residual(new Matrix1e(*res - *tmp));
      for (int i = 0; i != start->size(); ++i) residual->data(i) /= diag->data(i);
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      *prev -= *residual;
      if (::pow(residual->norm(),2.0) < tiny) break;
    }
    prev->print();
  }
#endif

  // testing BFGS update
  {
    cout << "  testing BFGS class" << endl;
    shared_ptr<Matrix1e> tmp(new Matrix1e(*target));

    // start with target/denom
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    BFGS<Matrix1e> bfgs(diag);

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 
      shared_ptr<Matrix1e> residual0(new Matrix1e(*res - *tmp));

      shared_ptr<Matrix1e> residual = bfgs.extrapolate(residual0, prev);

      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      *prev -= *residual;
      if (::pow(residual->norm(),2.0) < tiny) break;
    }
    prev->print();
  }

  // testing Linear with BFGS updates
  {
    cout << "  testing Linear class" << endl;
    shared_ptr<Matrix1e> tmp(new Matrix1e(*target));
    Linear<Matrix1e> linear(n*n, tmp);

    // start with target/denom
    shared_ptr<Matrix1e> prev(new Matrix1e(geom_));
    prev->element(0,0) = 1.0;

    BFGS<Matrix1e> bfgs(diag);

    for (int i = 0; i != n; ++i) {
      shared_ptr<Matrix1e> start(new Matrix1e(*prev));
      linear.orthog(start);
      shared_ptr<Matrix1e> res = start->clone();
      dgemv_("N", n*n, n*n, 1.0, hess.get(), n*n, start->data(), 1, 0.0, res->data(), 1); 

      shared_ptr<Matrix1e> residual0 = linear.compute_residual(start, res);

      shared_ptr<Matrix1e> residual = bfgs.extrapolate(residual0, linear.civec());
      cout << "residual " << setw(20) << setprecision(10) << fixed << residual->norm() << endl;
      if (::pow(residual->norm(),2.0) < tiny) break;

      prev = residual;
    }
    linear.civec()->print();
  }

#if 0
  // reference
  shared_ptr<Matrix1e> answer(new Matrix1e(*target));
  unique_ptr<int[]> ipiv(new int[n*n]);
  int info;
  dgesv_(n*n, 1, hess.get(), n*n, ipiv.get(), answer->data(), n*n, info);
  answer->print();
#endif
}
