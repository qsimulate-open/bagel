//
// BAGEL - Parallel electron correlation program.
// Filename: angularbatch.cc
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
// alxyz[1] later version.
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


#include <src/math/bessel.h>
#include <src/math/algo.h>
#include <src/integral/ecp/wigner3j.h>
#include <src/integral/ecp/angularbatch.h>
#include <src/integral/ecp/sphusplist.h>
#include <iomanip>

using namespace bagel;
using namespace std;

const static SphUSPList sphusplist;

double AngularBatch::integrate3SHs(array<pair<int, int>, 3> lm) const {

  const int l1 = lm[0].first;
  const int l2 = lm[1].first;
  const int l3 = lm[2].first;

  const static Wigner3j wigner3j;
  const double w1 = wigner3j.lookup_wigner3j(l1, lm[0].second, l2, lm[1].second, l3, lm[2].second);
  const double w2 = wigner3j.lookup_wigner3j(l1, 0, l2, 0, l3, 0);
  const double coeff = sqrt(0.25 * (2*l1 + 1) * (2*l2 + 1) * (2*l3 + 1) / pi__);

  return coeff * w1 * w2;

}

double AngularBatch::integrate3USP(array<int, 3> xyz_exponents) const {

  const int i = xyz_exponents[0];
  const int j = xyz_exponents[1];
  const int k = xyz_exponents[2];

  double out = 0.0;
  if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
    const static Factorial fact;
    const double num = 2.0 * fact(i) * fact(j) * fact(k) * fact((i+j+k+2)/2);
    const double denom = fact(i/2) * fact(j/2) * fact(k/2) * fact(i+j+k+2);
    out = 4.0 * pi__ * num / denom;
  }

  return out;

}

double AngularBatch::integrate2SH1USP(const pair<int, int> lm1, const pair<int, int> lm2, const array<int, 3> ijk) const {

  vector<double> usp1, usp2;
  sphusplist.sphuspfunc_call(lm1.first, lm1.second, usp1);
  sphusplist.sphuspfunc_call(lm2.first, lm2.second, usp2);

  double ans = 0.0;
  for (int i = 0; i != usp1.size(); ++i) {
    for (int j = 0; j != usp2.size(); ++j) {
      const double coeff = usp1[i] * usp2[j];
      if (coeff != 0.0) {
        array<int, 3> ki = {0, i, 0};
        for (int lp = lm1.first + 1; lp != 0; --lp)
          if (ki[1] - lp < 0) {
            ki[0] = lm1.first - ki[2] - ki[1];
            break;
          } else {
            ki[2]++;
            ki[1] -= lp;
          }
        array<int, 3> kj = {0, j, 0};
        for (int lp = lm2.first + 1; lp != 0; --lp)
          if (kj[1] - lp < 0) {
            kj[0] = lm1.first - kj[2] - kj[1];
            break;
          } else {
            kj[2]++;
            kj[1] -= lp;
          }
        const int x = ki[0] + kj[0] + ijk[0];
        const int y = ki[1] + kj[1] + ijk[1];
        const int z = ki[2] + kj[2] + ijk[2];
        array<int, 3> xyz = {x, y, z};
        ans += coeff * integrate3USP(xyz);
      }
    }
  }

  return ans;

}


double AngularBatch::project_one_centre(array<double, 3> posA, const array<int, 3> lxyz, const double expA,
                                        array<double, 3> posB, const array<int, 2> lm, const double r) {

  Comb comb;
  array<double, 3> AB;
  for (int i = 0; i != 3; ++i) AB[i] = posA[i] - posB[i];
  const double dAB = sqrt(AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]);
  const int nu = lxyz[0] + lxyz[1] + lxyz[2];
  const int lnu = lm[0] + nu;
  const double exponential = exp(-expA * (dAB - r) * (dAB - r)); // exp(-2aABr)i(2aABr) is computed, not i(2aABr)
  double ans = 0.0;
  for (int kx = 0; kx <= lxyz[0]; ++kx) {
    const double ckx = comb(lxyz[0], kx) * pow(AB[0], lxyz[0] - kx);
    for (int ky = 0; ky <= lxyz[1]; ++ky) {
      const double cky = comb(lxyz[1], ky) * pow(AB[1], lxyz[1] - ky);
      for (int kz = 0; kz <= lxyz[2]; ++kz) {
        const double ckz = comb(lxyz[2], kz) * pow(AB[2], lxyz[2] - kz);
        const int lk = kx + ky + kz;
        const double rxyz = pow(r, lk);

        double sld = 0.0;
        for (int ld = 0; ld <= lnu; ++ld) {
          double smu = 0.0;
          for (int m = 0; m <= 2 * ld; ++m) {
            const int mu = m - ld;
            shared_ptr<SphHarmonics> sphAB = make_shared<SphHarmonics>(ld, mu, AB);
            const double Z_AB = (dAB == 0 ? (1.0/sqrt(4.0*pi__)) : sphAB->zlm());

            const array<int, 3> exp = {kx, ky, kz};
            const pair<int, int> lm1(ld, mu);
            const pair<int, int> lm2(lm[0], lm[1]);
            smu += Z_AB * integrate2SH1USP(lm1, lm2, exp);
          }
          MSphBesselI msbessel(ld);
          const double sbessel = msbessel.compute(2.0 * expA * dAB * r);
          sld += smu * sbessel;
        }
        ans += sld * ckx * cky * ckz * rxyz * pow(-1.0, lk - nu);
      }
    }
  }

  return ans * 4.0 * pi__ * exponential;

}

double AngularBatch::project_many_centres(const double expA, const double expC, const double r) {

 vector<shared_ptr<const Shell_ECP>> shells_ecp = ecp_->shells_ecp();

 double ans = 0.0;
 for (auto& ishecp : shells_ecp) {
   const int l = ishecp->angular_number();
   if (l != ecp_->ecp_maxl()) {
     for (int i = 0; i != ishecp->ecp_exponents().size(); ++i) {
       if (ishecp->ecp_coefficients(i) != 0) {
         for (int m = 0; m <= 2*l; ++m) {
           array<int, 2> lm = {l, m - l};
           const double projA = project_one_centre(basisinfo_[0]->position(), ang0_, expA, ishecp->position(), lm, r);
           const double projC = project_one_centre(basisinfo_[1]->position(), ang1_, expC, ishecp->position(), lm, r);
           ans += ishecp->ecp_coefficients(i) * projA * projC *
                  pow(r, ishecp->ecp_r_power(i)) * exp(-ishecp->ecp_exponents(i) * r * r);
         }
       }
     }
   }
 }

 return ans;

}

double AngularBatch::compute(const double r) {

    double out = 0.0;

    const int begin0 = basisinfo_[0]->contraction_ranges(cont0_).first;
    const int end0   = basisinfo_[0]->contraction_ranges(cont0_).second;
    for (int i0 = begin0; i0 != end0; ++i0) {
      const double coef0 = basisinfo_[0]->contractions()[cont0_][i0];
      const double exp0  = basisinfo_[0]->exponents(i0);

      const int begin1 = basisinfo_[1]->contraction_ranges(cont1_).first;
      const int end1   = basisinfo_[1]->contraction_ranges(cont1_).second;
      for (int i1 = begin1; i1 != end1; ++i1) {
        const double coef1 = basisinfo_[1]->contractions()[cont1_][i1];
        const double exp1  = basisinfo_[1]->exponents(i1);

        out += coef0 * coef1 * project_many_centres(exp0, exp1, r);

      }
    }

  return out;

}

void AngularBatch::print() const {

  cout << "Compute the integral < shell_0 | lm > exp(-zeta r^n) < lm | shell_1> r^2 dr " << endl;
  cout << "Shell 0" << basisinfo_[0]->show() << endl;
  cout << "Shell 1" << basisinfo_[1]->show() << endl;
  cout << "ECP parameters" << endl;
  ecp_->print();

}

void AngularBatch::print_one_centre(array<double, 3> posA, const array<int, 3> lxyz, const double expA,
                                    array<double, 3> posB, const array<int, 2> lm, const double r) const {

  cout << "Project one centre <phiA | lmB> (r)" << endl;
  cout << "A = (" << posA[0] << ", " << posA[1] << ", " << posA[2] << ")" << endl;
  cout << "B = (" << posB[0] << ", " << posB[1] << ", " << posB[2] << ")" << endl;
  cout << "(kx, ky, kz) = (" << lxyz[0] << ", " << lxyz[1] << ", " << lxyz[2] << ")" << endl;
  cout << "(l, m)       = (" << lm[0] << ", " << lm[1] << ")" << endl;
  cout << "alpha = " << setw(15) << setprecision(9) << expA << endl;
  cout << "r     = " << setw(15) << setprecision(9) << r << endl;

}

