//
// BAGEL - Parallel electron correlation program.
// Filename: radial.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/integral/ecp/radial.h>

using namespace bagel;
using namespace std;

void RadialInt::integrate() {
  Timer radialtime;

  int n0 = 31;
  vector<int> sigma1(n0);
//transform_Becke(n0);
  transform_Ahlrichs(n0);
  vector<double> f = compute(r_);
  assert(f.size() == n0);
  for (int i = 0; i != r_.size(); ++i) sigma1[i] = i;

  double previous = inner_product(f.begin(), f.end(), w_.begin(), 0);
  int n1 = n0*2+1;

  for (int iter = 1; iter != max_iter_; ++iter) {
//  transform_Becke(n1); // Very slow!
//  transform_Log(n1, 3); //TODO: to be checked
    transform_Ahlrichs(n1);

    vector<int> sigma0(sigma1);
    sigma1.resize(n1);
    f.resize(n1);

    vector<double> rr(n0+1);
    for (int i = 0; i != n0; ++i) {
      sigma1[i] = sigma0[i]*2+1;
      sigma1[n0+i] = 2*i;
      rr[i] = r_[2*i];
    }
    sigma1[2*n0] = 2*n0;
    rr[n0] = r_[2*n0];

    vector<double> tmp = compute(rr);
    for (int i = 0; i <= n0; ++i) f[n0+i] = tmp[i];

    double ans = 0.0;
    for (int i = 0; i != r_.size(); ++i) ans += f[i] * w_[sigma1[i]];

    const double error = ans - previous;
    if (print_intermediate_)
       cout << "Iter = " << setw(5) << iter << setw(10) << "npts = " << setw(10) << n1
                 << setw(10) << "ans = " << setw(20) << setprecision(10) << ans
                 << setw(10) << "err = " << setw(20) << setprecision(10) << error << endl;
    if (fabs(error) < thresh_int_) {
      if (print_intermediate_) {
        cout << "Integration converged..." << endl;
        cout << "Radial integral = " << ans << endl;
        cout << "Radial time = " << radialtime.tick() << endl;
      }
      integral_ = ans;
      break;
    } else if (iter == max_iter_-1) {
      cout << "Max iteration exceeded..." << endl;
    }
    previous = ans;
    x_.clear();
    w_.clear();
    r_.clear();
    n0 = n1;
    n1 = 2*n1+1;
  }
}

void RadialInt::transform_Log(const int ngrid, const int m) { // Mura and Knowles JCP, 104, 9848.
  w_.resize(ngrid);
  r_.resize(ngrid);
  const double alpha = 5.0;
  for (int i = 1; i <= ngrid; ++i) {
    const double x = i / (ngrid + 1.0);
    const double xm = 1.0 - pow(x, m);
    r_[i-1] = - alpha * log(xm);
    w_[i-1] = pow(r_[i-1], 2) * alpha * m * pow(x, m-1) / (xm * (ngrid + 1.0));
  }
}

void RadialInt::transform_Ahlrichs(const int ngrid) { // Treutler and Ahlrichs JCP, 102, 346.
  GaussChebyshev2nd(ngrid);
  r_.resize(ngrid);
  const double alpha = 1.0;
  for (int i = 0; i != ngrid; ++i) {
    const double exp = 0.6;
    const double prefactor = alpha / log(2.0);
    r_[i]  = prefactor * pow(1.0 + x_[i], exp) * log(2.0 / (1 - x_[i]));
    w_[i] *= prefactor * (exp * pow(1.0 + x_[i], exp - 1.0) * log(2.0 / (1.0 - x_[i]))
                      + pow(1.0 + x_[i], exp) / (1.0 - x_[i]));
  }
}

void RadialInt::transform_Becke(const int ngrid) { // Becke JCP, 88, 2547.
  GaussChebyshev2nd(ngrid);
  r_.resize(ngrid);
  const double alpha = 1.0;
  for (int i = 0; i != ngrid; ++i) {
    r_[i] = alpha * (1.0 + x_[i]) / (1.0 - x_[i]);
    w_[i] *= 2.0 * alpha / pow(1.0 - x_[i], 2);
  }
}

void RadialInt::GaussChebyshev2nd(const int ngrid) {
  x_.resize(ngrid);
  w_.resize(ngrid);
  for (int i = 1; i <= ngrid; ++i) {
    x_[i-1] = cos(i * pi__ / (ngrid + 1));
    w_[i-1] = pi__ * sin(i * pi__ / (ngrid + 1)) / (ngrid + 1);
  }
}
