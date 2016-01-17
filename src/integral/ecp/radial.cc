//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: radial.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/integral/ecp/radial.h>

using namespace bagel;
using namespace std;

void RadialInt::integrate() {
  Timer radialtime;

  integral_.resize(nc_);

  int n0 = 31;
  vector<int> sigma1(n0);
//transform_Becke(n0);
  transform_Ahlrichs(n0);
  assert(r_.size() == n0);
  vector<double> f = compute(r_);
  assert(f.size() == n0*nc_);
  for (int i = 0; i != n0; ++i) sigma1[i] = i;

  vector<double> previous(nc_, 0.0);
  for (int ic = 0; ic != nc_; ++ic)
    for (int i = 0; i != n0; ++i) previous[ic] = f[ic*n0+i] * w_[i];
    //previous[ic] = inner_product(f.begin()+ic*n0, f.begin()+(ic+1)*n0-1, w_.begin(), 0);

  int n1 = n0*2+1;

  for (int iter = 1; iter != max_iter_; ++iter) {
//  transform_Becke(n1); // Very slow!
//  transform_Log(n1, 3); //TODO: to be checked
    transform_Ahlrichs(n1);

    vector<int> sigma0(sigma1);
    sigma1.resize(n1);
    f.resize(n1*nc_);
    for (int ic = nc_-1; ic > 0; --ic)
      for (int i = 0; i < n0; ++i) f[ic*n1+i] = f[ic*n0+i];

    vector<double> rr(n0+1);
    for (int i = 0; i != n0; ++i) {
      sigma1[i] = sigma0[i]*2+1;
      sigma1[n0+i] = 2*i;
      rr[i] = r_[2*i];
    }
    sigma1[2*n0] = 2*n0;
    rr[n0] = r_[2*n0];

    vector<double> tmp = compute(rr);
    for (int ic = 0; ic != nc_; ++ic)
      for (int i = 0; i <= n0; ++i) f[ic*n1+n0+i] = tmp[ic*(n0+1)+i];

    vector<double> ans(nc_, 0.0);
    for (int ic = 0; ic != nc_; ++ic)
      for (int i = 0; i != n1; ++i) ans[ic] += f[ic*n1+i] * w_[sigma1[i]];

    vector<double> error(nc_, 0.0);
    double maxerror = 0.0;
    for (int ic = 0; ic != nc_; ++ic) {
      error[ic] = fabs(ans[ic] - previous[ic]);
      if (error[ic] > maxerror) maxerror = error[ic];
    }
    if (print_intermediate_) {
       cout << "Iter = " << setw(5) << iter << setw(10) << "npts = " << setw(10) << n1 << endl;
       for (int ic = 0; ic != nc_; ++ic)
         cout << setw(10) << "ans[" << ic << "] = " << setw(20) << setprecision(12) << ans[ic]
              << setw(10) << "err[" << ic << "] = " << setw(20) << setprecision(12) << error[ic] << endl;
    }
    if (maxerror <= thresh_int_) {
      if (print_intermediate_) {
        cout << "Integration converged in " << iter << " interations:   ";
        for (int ic = 0; ic != nc_; ++ic) {
          cout << setw(20) << setprecision(12) << ans[ic];
          if (nc_ > 1) cout << ",  ";
        }
        cout << endl;
        radialtime.tick_print("ECP radial integration");
      }
      for (int ic = 0; ic != nc_; ++ic)  integral_[ic] = ans[ic];
      break;
    } else if (iter == max_iter_-1) {
      throw runtime_error("Max iteration exceeded in ecp/radial...");
    }
    for (int ic = 0; ic != nc_; ++ic) previous[ic] = ans[ic];
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
