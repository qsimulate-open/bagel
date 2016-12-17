//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: soangularbatch.cc
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
#include <src/integral/ecp/soangularbatch.h>
#include <src/integral/ecp/sphusplist.h>
#include <iomanip>

using namespace bagel;
using namespace std;

const static SphUSPList sphusplist;
const static DoubleFactorial df;

SOBatch::SOBatch(const shared_ptr<const SOECP> _so, const array<shared_ptr<const Shell>,2>& _info,
                           const int contA, const int contC, const array<int, 3> angA, const array<int, 3> angC,
                           const bool print, const int max_iter, const double thresh_int)
 : RadialInt(3, print, max_iter, thresh_int),
   basisinfo_(_info), so_(_so), cont0_(contA), cont1_(contC), ang0_(angA), ang1_(angC) {

  init();
  map_angular_number();

}

void SOBatch::map_angular_number() {

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

complex<double> SOBatch::theta(const int m) const {
  const double tau = (m < 0) ? 0.0 : 1.0;

  const double re = (m==0) ? 0.5 : tau/sqrt(2.0);
  const double im = (m==0) ? 0.0 : (1.0-tau)/sqrt(2.0);

  return complex<double>(re, im);
}

array<double, 3> SOBatch::fm0lm1(const int l, const int m0, const int m1) const { /* Im{aa}, Re{ab}, Im{ab} */

  assert(l>0 && abs(m0)<=l && abs(m1)<=l && m1<m0);
  array<double, 3> out = {{0.0, 0.0, 0.0}};

  if (m0 == -m1) out[0] = 0.5 * m1;

  const double deltap = (abs(m0) == abs(m1) + 1) ? 1.0 : 0.0;
  const double deltam = (abs(m0) == abs(m1) - 1) ? 1.0 : 0.0;

  const double tau0 = (m0 < 0) ? 0.0 : 1.0;
  const double tau1 = (m1 < 0) ? 0.0 : 1.0;

  complex<double> mu = conj(theta(m0))*theta(m1);
  if (m0*m1==0) mu *= 2.0;

  complex<double> ab = 0.5 * mu * (deltap - pow(-1.0, tau0+tau1) * deltam);
  ab *= sqrt((l + m0*m0 - abs(m0*m1)) * (l + m1*m1 - abs(m0*m1)));

  out[1] = real(ab);
  out[2] = imag(ab);

  return out;
}

double SOBatch::angularA(const int h, const int ld, const vector<double> usp) {

  double out = 0;

  const int l = static_cast<int>(round(sqrt(usp.size()*2)-1));
  assert((l+1)*(l+2) == usp.size() * 2);

  for (int a = max(0, h-ang0_[1]-ang0_[2]); a <= min(ang0_[0], h); ++a) {
    for (int b = max(0, h-a-ang0_[2]); b <= min(ang0_[1], h-a); ++b) {
      const int index = a * ANG_HRR_END * ANG_HRR_END + b * ANG_HRR_END + h - a - b;
      if (abs(c0_[index]) > 1e-15) {

        double smu = 0.0;
        for (int j = 0; j != usp.size(); ++j) {
          if (usp[j] != 0.0) {
            map<int, array<int, 3>>::const_iterator pj = map_[l].find(j);
            assert (pj != map_[l].end());
            const array<int, 3> kj = pj->second;

            for (int mu = 0; mu <= 2*ld; ++mu) {
              const vector<double> usp1 = sphusplist.sphuspfunc_call(ld, mu-ld);
              double sAB = 0.0;
              for (int i = 0; i != usp1.size(); ++i) {
                if (usp1[i] != 0.0) {
                  map<int, array<int, 3>>::const_iterator pi = map_[ld].find(i);
                  assert (pi != map_[ld].end());
                  const array<int, 3> ki = pi->second;
                  const int x = ki[0] + kj[0] + a;
                  const int y = ki[1] + kj[1] + b;
                  const int z = ki[2] + kj[2] + h-a-b;
                  const double xyz = (x % 2 == 0 && y % 2 == 0 && z % 2 == 0) ? (4.0 * pi__ * df(x-1) * df(y-1) * df(z-1) / df(x+y+z+1)) : 0.0;
                  sAB += usp1[i] * usp[j] * xyz;
                }
              }
              smu += zAB_[ld][mu] * sAB;
            }
          }
        }
        out += smu*c0_[index];
      }
    }
  }

  return out;

}

double SOBatch::angularC(const int h, const int ld, const vector<double> usp) {

  double out = 0;

  const int l = static_cast<int>(round(sqrt(usp.size()*2)-1));
  assert((l+1)*(l+2) == usp.size() * 2);

  for (int a = max(0, h-ang1_[1]-ang1_[2]); a <= min(ang1_[0], h); ++a) {
    for (int b = max(0, h-a-ang1_[2]); b <= min(ang1_[1], h-a); ++b) {
      const int index = a * ANG_HRR_END * ANG_HRR_END + b * ANG_HRR_END + h - a - b;
      if (abs(c1_[index]) > 1e-15) {

        double smu = 0.0;
        for (int j = 0; j != usp.size(); ++j) {
          if (usp[j] != 0.0) {
            map<int, array<int, 3>>::const_iterator pj = map_[l].find(j);
            assert (pj != map_[l].end());
            const array<int, 3> kj = pj->second;

            for (int mu = 0; mu <= 2*ld; ++mu) {
              const vector<double> usp1 = sphusplist.sphuspfunc_call(ld, mu-ld);
              double sCB = 0.0;
              for (int i = 0; i != usp1.size(); ++i) {
                if (usp1[i] != 0.0) {
                  map<int, array<int, 3>>::const_iterator pi = map_[ld].find(i);
                  assert (pi != map_[ld].end());
                  const array<int, 3> ki = pi->second;
                  const int x = ki[0] + kj[0] + a;
                  const int y = ki[1] + kj[1] + b;
                  const int z = ki[2] + kj[2] + h-a-b;
                  const double xyz = (x % 2 == 0 && y % 2 == 0 && z % 2 == 0) ? (4.0 * pi__ * df(x-1) * df(y-1) * df(z-1) / df(x+y+z+1)) : 0.0;
                  sCB += usp1[i] * usp[j] * xyz;
                }
              }
              smu += zCB_[ld][mu] * sCB;
            }
          }
        }
        out += smu*c1_[index];
      }
    }
  }

  return out;

}

vector<double> SOBatch::project(const int l, const vector<double> r) {

  const static MSphBesselI msbessel;
  vector<vector<double>> rbessel0(r.size()), rbessel1(r.size());

  const int begin0 = basisinfo_[0]->contraction_ranges(cont0_).first;
  const int end0   = basisinfo_[0]->contraction_ranges(cont0_).second;

  const int begin1 = basisinfo_[1]->contraction_ranges(cont1_).first;
  const int end1   = basisinfo_[1]->contraction_ranges(cont1_).second;

  for (int ir = 0; ir != r.size(); ++ir) {
    vector<double> b0(l0_+l+1, 0.0);
    for (int i0 = begin0; i0 != end0; ++i0) {
      const double coef0 = basisinfo_[0]->contractions()[cont0_][i0];
      const double exp0 = basisinfo_[0]->exponents(i0);
      const double fac0 = coef0 * exp(-exp0 * pow(dAB_-r[ir], 2));
      for (int i = 0; i <= l0_+l; ++i)
        b0[i] += fac0 * msbessel.compute(i, 2.0*exp0*dAB_*r[ir]);
    }
    rbessel0[ir] = b0;

    vector<double> b1(l1_+l+1, 0.0);
    for (int i1 = begin1; i1 != end1; ++i1) {
      const double coef1 = basisinfo_[1]->contractions()[cont1_][i1];
      const double exp1 = basisinfo_[1]->exponents(i1);
      const double fac1 = coef1 * exp(-exp1 * pow(dCB_-r[ir], 2));
      for (int i = 0; i <= l1_+l; ++i)
        b1[i] += fac1 * msbessel.compute(i, 2.0*exp1*dCB_*r[ir]);
    }
    rbessel1[ir] = b1;
  }

  vector<vector<double>> usp(2*l+1);
  for (int m = 0; m <= 2*l; ++m) usp[m] = sphusplist.sphuspfunc_call(l, m-l);

  vector<double> out(3*r.size(), 0.0);

  for (int ld0 = max(0, l-l0_); ld0 <= l+l0_; ++ld0) {
    for (int ld1 = max(0, l-l1_); ld1 <= l+l1_; ++ld1) {
      const int c0 = l0_ - (l0_ - abs(ld0-l))%2;
      const int c1 = l1_ - (l1_ - abs(ld1-l))%2;

      const int gmin = abs(ld0-l)+abs(ld1-l);
      const int gmax = c0 + c1;
      for (int g = gmin; g <= gmax; g += 2) {
        array<double, 3> sum = {{0.0, 0.0, 0.0}};

        const int hmin = max(abs(ld0-l), g - c1);
        const int hmax = min(c0, g - abs(ld1-l));
        for (int i = 0; i != fm0lm1_[l-1].size(); ++i) {
          tuple<int, int, int, double> fmm = fm0lm1_[l-1][i];
          const int id = get<0>(fmm);
          const int m0 = get<1>(fmm);
          const int m1 = get<2>(fmm);
          const double f = get<3>(fmm);
          for (int h = hmin; h <= hmax; h += 2)
            sum[id] += (angularA(h, ld0, usp[m0]) * angularC(g-h, ld1, usp[m1]) - angularA(h, ld0, usp[m1]) * angularC(g-h, ld1, usp[m0]))*f;
        }
        for (int ir = 0; ir != r.size(); ++ir) {
          const double p = rbessel0[ir][ld0] * rbessel1[ir][ld1] * pow(r[ir], g);
          for (int id = 0; id != 3; ++id) {
            const int index = id*r.size() + ir;
            out[index] += sum[id] * p;
          }
        }

      }
    }
  }

  return out;

}

vector<double> SOBatch::compute(const vector<double> r) {

  vector<shared_ptr<const Shell_ECP>> shells_so = so_->shells_so();
  vector<double> out(3*r.size(), 0.0);

  for (auto& ishso : shells_so) {
    const int l = ishso->angular_number();
    vector<double> p = project(l, r);
    for (int i = 0; i != ishso->ecp_exponents().size(); ++i)
      if (ishso->ecp_coefficients(i) != 0) {
        const double ecpcoeff = 16.0 * pi__ * pi__ * ishso->ecp_coefficients(i); // 2/(2l+1) not needed for Stuggart basis sets
        for (int ir = 0; ir != r.size(); ++ir) {
          const double coeff = ecpcoeff * pow(r[ir], ishso->ecp_r_power(i)) * exp(-ishso->ecp_exponents(i) * r[ir] * r[ir]);
          for (int id = 0; id != 3; ++id) {
            const int index = id*r.size() + ir;
            out[index] += coeff * p[index];
          }
        }
      }
  }

  return out;

}

void SOBatch::init() {

  l0_ = ang0_[0] + ang0_[1] + ang0_[2];
  l1_ = ang1_[0] + ang1_[1] + ang1_[2];

  for (int i = 0; i != 3; ++i) {
    AB_[i] = basisinfo_[0]->position(i) - so_->position(i);
    CB_[i] = basisinfo_[1]->position(i) - so_->position(i);
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
        c0_[index] = ckx * cky * ckz * pow(-1.0, l0_-kx-ky-kz);
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
        c1_[index] = ckx * cky * ckz * pow(-1.0, l1_-kx-ky-kz);
      }
    }
  }

  for (int l = 0; l <= max(l0_, l1_) + so_->so_maxl(); ++l) {
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

  fm0lm1_.resize(so_->so_maxl());
  for (int l = 1; l <= so_->so_maxl(); ++l) {
    vector<tuple<int, int, int, double>> fmm;
    fmm.resize((2*l-1)*3);
    int index = 0;
    for (int m0 = 0; m0 <= 2*l; ++m0)
    for (int m1 = 0; m1 <= m0-1; ++m1) {
      const array<double, 3> f = fm0lm1(l, m0-l, m1-l);
      for (int i = 0; i < 3; ++i) if (abs(f[i]) > 1e-15) {
        fmm[index] = make_tuple(i, m0, m1, f[i]);
        ++index;
      }
    }
    fmm.resize(index);
    fm0lm1_[l-1] = fmm;
    assert (index <= (2*l-1)*3);
  }

}

void SOBatch::print() const {
  cout << "+++ Shell 0 +++" << endl << basisinfo_[0]->show() << endl;
  cout << "+++ Shell 1 +++" << endl << basisinfo_[1]->show() << endl;
  so_->print();
}
