// Author: Hai-Anh Le
// Date: Jul 30, 2014

#include <map>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <src/util/math/sphharmonics.h>
#include "src/integral/ecp/sphusplist.h"
#include "cartgauss.h"

using namespace std;
using namespace bagel;
using namespace test;

const static DoubleFactorial df;
const static SphUSPList list;

namespace test {

class SOInt {
  protected:
    shared_ptr<const CartesianGauss> carA_;
    shared_ptr<const SphHarmonics> sphB_;
    array<int, 3> ang0_;
    int l0_;
    array<double, 3> AB_;
    double dAB_;
    vector<double> c0_;
    vector<std::vector<double>> zAB_;
    vector<std::map<int, std::array<int, 3>>> map_;

  public:
    SOInt(shared_ptr<const CartesianGauss> c, shared_ptr<const SphHarmonics> sh)
      : carA_(c), sphB_(sh) { init(); }

    ~SOInt() {}

    void init() {
      ang0_ = carA_->angular_momentum();
      l0_ = ang0_[0] + ang0_[1] + ang0_[2];
      for (int i = 0; i != 3; ++i) {
        AB_[i] = carA_->centre(i) - sphB_->centre(i);
      }
      dAB_ = sqrt(AB_[0]*AB_[0] + AB_[1]*AB_[1] + AB_[2]*AB_[2]);

      map_angular_number();
      compute_c0();
      compute_zAB();
    }

    void print() const {
      carA_->print();
      sphB_->print();
    }

    void map_angular_number() {
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

    double angularA(const int h, const int ld, const int l, const vector<double> usp) {
      double out = 0;

      //const int l = static_cast<int>(round(sqrt(usp.size()*2)-1));
      //assert((l+1)*(l+2) == static_cast<int>(usp.size() * 2));

      const int jend = (l+1)*(l+2)/2;
      const int iend = (ld+1)*(ld+2)/2;

      const int amin = max(0, h-ang0_[1]-ang0_[2]);
      const int amax = min(ang0_[0], h);
      for (int a = amin; a <= amax; ++a) {
        const int bmin = max(0, h-a-ang0_[2]);
        const int bmax = min(ang0_[1], h-a);
        for (int b = bmin; b <= bmax; ++b) {
          const int index = a * ANG_HRR_END * ANG_HRR_END + b * ANG_HRR_END + h - a - b;
          if (c0_[index] > 1e-15) {
            double smu = 0.0;
            for (int j = 0; j != jend; ++j) {
              if (usp[j] != 0.0) {
                map<int, array<int, 3>>::const_iterator pj = map_[l].find(j);
                assert (pj != map_[l].end());
                const array<int, 3> kj = pj->second;

                for (int mu = 0; mu <= 2*ld; ++mu) {
                  const vector<double> usp1 = list.sphuspfunc_call(ld, mu-ld);
                  double sAB = 0.0;
                  for (int i = 0; i != iend; ++i) {
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

    void compute_c0() {

      c0_.resize(ANG_HRR_END*ANG_HRR_END*ANG_HRR_END);

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
    }

    void compute_zAB() {
      for (int l = 0; l <= l0_ + sphB_->angular_momentum(0); ++l) {
        vector<double> zAB_l(2*l+1, 0.0);
        for (int m = 0; m <= 2*l; ++m) {
          auto shAB = make_shared<SphHarmonics>(l, m-l, AB_);
          zAB_l[m] = (dAB_ < 1e-12 ? (1.0/sqrt(4.0*pi__)) : shAB->zlm());
        }
        zAB_.push_back(zAB_l);
      }
    }
};

}

using namespace test;
int main() {

  cout << " +++ TEST P(h, ld, l, m) +++ " << endl;

  const array<double, 3> centreA = {{0.0000, 0.0000, 0.0000}};
  const array<double, 3> centreB = {{0.0000, 0.0000, 0.0000}};
  const double alphaA = 1.0;

  const int jAmax = 1;
  for (int jA = 1; jA <= jAmax; ++jA) {
    for (int lz = 0; lz <= jA; ++lz) {
      for (int ly = 0; ly <= jA-lz; ++ly) {
        const int lx = jA-lz-ly;
        const array<int, 3> angA = {{lx, ly, lz}};
        auto cargaussA = make_shared<const CartesianGauss>(alphaA, angA, centreA);

        const int lmax = 1;
        for (int l = 0; l <= lmax; ++l) {
          for (int m = 0; m <= 2*l; ++m) {
            const array<int, 2> lm = {{l, m-l}};
            auto rsh = make_shared<const SphHarmonics>(lm, centreB);
            const vector<double> usp = list.sphuspfunc_call(l, m-l);
            auto so = make_shared<SOInt>(cargaussA, rsh);

            for (int ld = 0; ld <= l+jA; ++ld) {
              for (int h = 0; h <= jA; ++h) {
                const double p = so->angularA(h, ld, l, usp);
                if (p > 1e-12) {
                  cout << "h = " << h << " (lx, ly, lz) = (" << lx << ", " << ly << ", " << lz << ")" << endl;
                  cout << "P(" << h << ", " << ld << ", " << l << ", " << m-l << ") = " << setprecision(12) << p << endl;
                  cout << endl;
                }
              }
            }
          }
        }
      }
    }
  }

  return 0;
}
