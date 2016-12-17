//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: angularbatch.cc
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


#include <src/util/math/bessel.h>
#include <src/util/math/algo.h>
#include <src/integral/ecp/wigner3j.h>
#include <src/integral/ecp/angularbatch.h>
#include <src/integral/ecp/sphusplist.h>
#include <iomanip>

using namespace bagel;
using namespace std;

const static SphUSPList sphusplist;
const static DoubleFactorial df;

AngularBatch::AngularBatch(const shared_ptr<const ECP> _ecp, const array<shared_ptr<const Shell>,2>& _info,
                           const int contA, const int contC, const array<int, 3> angA, const array<int, 3> angC,
                           const bool print, const int max_iter, const double thresh_int)
 : RadialInt(1, print, max_iter, thresh_int),
   basisinfo_(_info), ecp_(_ecp), cont0_(contA), cont1_(contC), ang0_(angA), ang1_(angC) {

  init();
  map_angular_number();

}

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

void AngularBatch::map_angular_number() {

  for (int l = 0; l != ANG_HRR_END*2-1; ++l) {
    map<int, array<int, 3>> mapl;
    int key = 0;
    for (int z = 0; z <= l; ++z) {
      for (int y = 0; y <= l - z; ++y) {
        const int x = l - y - z;
        array<int, 3> xyz = {{x, y, z}};
        mapl.insert(make_pair(key, xyz));
        ++key;
      }
    }
    map_.push_back(mapl);
  }

}

vector<double> AngularBatch::project_AB(const int l, const vector<double> usp, const vector<double> r) {
  const static MSphBesselI msbessel;

  vector<vector<double>> rbessel(r.size());

  const int begin0 = basisinfo_[0]->contraction_ranges(cont0_).first;
  const int end0   = basisinfo_[0]->contraction_ranges(cont0_).second;
  for (int ir = 0; ir != r.size(); ++ir) {
    vector<double> bessel(l0_+l+1);
    for (int i0 = begin0; i0 != end0; ++i0) {
      const double coef0 = basisinfo_[0]->contractions()[cont0_][i0];
      const double exp0  = basisinfo_[0]->exponents(i0);
      const double fac = coef0 * exp(-exp0 * pow(dAB_-r[ir], 2));
      for (int i = 0; i <= l0_+l; ++i)
        bessel[i] += fac * msbessel.compute(i, 2.0 * exp0 * dAB_ * r[ir]);
    }
    rbessel[ir] = bessel;
  }

  vector<double> out(r.size(), 0.0);
  for (int j = 0; j != usp.size(); ++j) {
    if (usp[j] != 0.0) {
      map<int, array<int, 3>>::const_iterator pj = map_[l].find(j);
      assert (pj != map_[l].end());
      const array<int, 3> kj = pj->second;

      for (int kx = 0; kx <= ang0_[0]; ++kx)
      for (int ky = 0; ky <= ang0_[1]; ++ky)
      for (int kz = 0; kz <= ang0_[2]; ++kz) {
        const int lk = kx + ky + kz;
        const int index = kx * ANG_HRR_END * ANG_HRR_END + ky * ANG_HRR_END + kz;
        const double coeff = c0_[index] * pow(-1.0, lk - l0_);
        if (abs(coeff) > 1e-15) {
        for (int ld = max(l-lk, 0); ld <= l+lk; ++ld) {
            if ((l + lk - ld) % 2 == 0) {
              double smu = 0.0;
              for (int mu = 0; mu <= 2 * ld; ++mu) {

                const vector<double> usp1 = sphusplist.sphuspfunc_call(ld, mu-ld);
                double sAB = 0.0;
                for (int i = 0; i != usp1.size(); ++i) {
                  if (usp1[i] != 0.0) {
                    map<int, array<int, 3>>::const_iterator p = map_[ld].find(i);
                    assert (p != map_[ld].end());
                    const array<int, 3> ki = p->second;
                    const int x = ki[0] + kj[0] + kx;
                    const int y = ki[1] + kj[1] + ky;
                    const int z = ki[2] + kj[2] + kz;
                    const double xyz = (x % 2 == 0 && y % 2 == 0 && z % 2 == 0) ? (4.0 * pi__ * df(x-1) * df(y-1) * df(z-1) / df(x+y+z+1)) : 0.0;
                    sAB += usp1[i] * usp[j] * xyz;
                  }
                }
                smu += zAB_[ld][mu] * sAB;
              }
              for (int ir = 0; ir != r.size(); ++ir) out[ir] += smu * rbessel[ir][ld] * coeff * pow(r[ir], lk);
            }
          }
        }
      }
    }
  }

  return out;

}

vector<double> AngularBatch::project_CB(const int l, const vector<double> usp, const vector<double> r) {
  const static MSphBesselI msbessel;

  vector<vector<double>> rbessel(r.size());

  const int begin1 = basisinfo_[1]->contraction_ranges(cont1_).first;
  const int end1   = basisinfo_[1]->contraction_ranges(cont1_).second;
  for (int ir = 0; ir != r.size(); ++ir) {
    vector<double> bessel(l1_+l+1);
    for (int i1 = begin1; i1 != end1; ++i1) {
      const double coef1 = basisinfo_[1]->contractions()[cont1_][i1];
      const double exp1  = basisinfo_[1]->exponents(i1);
      const double fac = coef1 * exp(-exp1 * pow(dCB_-r[ir], 2));
      for (int i = 0; i <= l1_+l; ++i)
        bessel[i] += fac * msbessel.compute(i, 2.0 * exp1 * dCB_ * r[ir]);
    }
    rbessel[ir] = bessel;
  }

  vector<double> out(r.size(), 0.0);
  for (int j = 0; j != usp.size(); ++j) {
    if (usp[j] != 0.0) {
      map<int, array<int, 3>>::const_iterator pj = map_[l].find(j);
      assert (pj != map_[l].end());
      const array<int, 3> kj = pj->second;

      for (int kx = 0; kx <= ang1_[0]; ++kx)
      for (int ky = 0; ky <= ang1_[1]; ++ky)
      for (int kz = 0; kz <= ang1_[2]; ++kz) {
        const int lk = kx + ky + kz;
        const int index = kx * ANG_HRR_END * ANG_HRR_END + ky * ANG_HRR_END + kz;
        const double coeff = c1_[index] * pow(-1.0, lk - l1_);
        if (abs(coeff) > 1e-15) {
          for (int ld = max(l-lk, 0); ld <= l+lk; ++ld) {
            if ((l + lk - ld) % 2 == 0) {
              double smu = 0.0;
              for (int mu = 0; mu <= 2 * ld; ++mu) {

                const vector<double> usp1 = sphusplist.sphuspfunc_call(ld, mu-ld);
                double sCB = 0.0;
                for (int i = 0; i != usp1.size(); ++i) {
                  if (usp1[i] != 0.0) {
                    map<int, array<int, 3>>::const_iterator p = map_[ld].find(i);
                    assert (p != map_[ld].end());
                    const array<int, 3> ki = p->second;
                    const int x = ki[0] + kj[0] + kx;
                    const int y = ki[1] + kj[1] + ky;
                    const int z = ki[2] + kj[2] + kz;
                    const double xyz = (x % 2 == 0 && y % 2 == 0 && z % 2 == 0) ? (4.0 * pi__ * df(x-1) * df(y-1) * df(z-1) / df(x+y+z+1)) : 0.0;
                    sCB += usp1[i] * usp[j] * xyz;
                  }
                }
                smu += zCB_[ld][mu] * sCB;
              }
              for (int ir = 0; ir != r.size(); ++ir) out[ir] += smu * rbessel[ir][ld] * coeff * pow(r[ir], lk);
            }
          }
        }
      }
    }
  }

  return out;

}

vector<double> AngularBatch::compute(const vector<double> r) {

  vector<shared_ptr<const Shell_ECP>> shells_ecp = ecp_->shells_ecp();
  vector<double> out(r.size(), 0.0);

  for (auto& ishecp : shells_ecp) {
    const int l = ishecp->angular_number();
    if (l != ecp_->ecp_maxl()) {
      for (int m = 0; m <= 2*l; ++m) {
        const vector<double> usp = sphusplist.sphuspfunc_call(l, m-l);
        vector<double> pA = project_AB(l, usp, r);
        vector<double> pC = project_CB(l, usp, r);
        for (int i = 0; i != ishecp->ecp_exponents().size(); ++i)
          if (ishecp->ecp_coefficients(i) != 0) {
            const double coeff = 16.0 * pi__ * pi__ * ishecp->ecp_coefficients(i);
            for (int ir = 0; ir != r.size(); ++ir)
              out[ir] += coeff * pow(r[ir], ishecp->ecp_r_power(i)) * exp(-ishecp->ecp_exponents(i) * r[ir] * r[ir]) * pA[ir] * pC[ir];
          }
      }
    }
  }


  return out;

}

void AngularBatch::init() {

  for (int i = 0; i != 3; ++i) {
    AB_[i] = basisinfo_[0]->position(i) - ecp_->position(i);
    CB_[i] = basisinfo_[1]->position(i) - ecp_->position(i);
  }
  dAB_ = sqrt(pow(AB_[0], 2) + pow(AB_[1], 2) + pow(AB_[2], 2));
  dCB_ = sqrt(pow(CB_[0], 2) + pow(CB_[1], 2) + pow(CB_[2], 2));

  c0_.resize(ANG_HRR_END*ANG_HRR_END*ANG_HRR_END);
  c1_.resize(ANG_HRR_END*ANG_HRR_END*ANG_HRR_END);

  const static Comb c;

  for (int kx = 0; kx <= ang0_[0]; ++kx) {
    const double ckx = c(ang0_[0], kx) * pow(AB_[0], ang0_[0] - kx);
    for (int ky = 0; ky <= ang0_[1]; ++ky) {
      const double cky = c(ang0_[1], ky) * pow(AB_[1], ang0_[1] - ky);
      for (int kz = 0; kz <= ang0_[2]; ++kz) {
        const double ckz = c(ang0_[2], kz) * pow(AB_[2], ang0_[2] - kz);
        const int index = kx * ANG_HRR_END * ANG_HRR_END + ky * ANG_HRR_END + kz;
        c0_[index] = ckx * cky * ckz;
      }
    }
  }

  for (int kx = 0; kx <= ang1_[0]; ++kx) {
    const double ckx = c(ang1_[0], kx) * pow(CB_[0], ang1_[0] - kx);
    for (int ky = 0; ky <= ang1_[1]; ++ky) {
      const double cky = c(ang1_[1], ky) * pow(CB_[1], ang1_[1] - ky);
      for (int kz = 0; kz <= ang1_[2]; ++kz) {
        const double ckz = c(ang1_[2], kz) * pow(CB_[2], ang1_[2] - kz);
        const int index = kx * ANG_HRR_END * ANG_HRR_END + ky * ANG_HRR_END + kz;
        c1_[index] = ckx * cky * ckz;
      }
    }
  }

  l0_ = ang0_[0] + ang0_[1] + ang0_[2];
  l1_ = ang1_[0] + ang1_[1] + ang1_[2];

  for (int l = 0; l != max(l0_, l1_) + ecp_->ecp_maxl(); ++l) {
    vector<double> zAB_l(2*l+1, 0.0), zCB_l(2*l+1, 0.0);
    for (int m = 0; m <= 2*l; ++m) {
      auto shAB = make_shared<SphHarmonics>(l, m-l, AB_);
      zAB_l[m] = (dAB_ < 1e-12 ? (1.0/sqrt(4.0*pi__)) : shAB->zlm());
      auto shCB = make_shared<SphHarmonics>(l, m-l, CB_);
      zCB_l[m] = (dCB_ < 1e-12 ? (1.0/sqrt(4.0*pi__)) : shCB->zlm());
    }
    zAB_.push_back(zAB_l);
    zCB_.push_back(zCB_l);
  }

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

