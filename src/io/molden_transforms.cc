//
// Author : Toru Shiozaki
// Date   : May 2009
//
// Modified by :  Shane Parker
// Date        :  July 2012
//
// This was originally code found in src/rysint/carsph_gen/carsph_gen.cc
// but repurposed here for the molden section. At the moment, it is limited to at most
// f functions.
//

#include <src/io/moldenin.h>

using namespace std;
using namespace bagel;

#define LARGE 32
#define LEND 5


void MoldenIn::compute_transforms() {
  struct Data {
    vector<double> factorial;
    Data() {
      factorial.push_back(1.0);
      for (int i = 1; i != LEND * 2; ++i)
        factorial.push_back(factorial.back() * i);
    }
    double comb(const int i, const int j) const {
      return factorial[i] / factorial[j] / factorial[i - j];
    }
  } data;

  const double one = 1.0;
  const double two = 2.0;
  const double quarter = 0.25;

  vector<double> factorial = data.factorial;

  vector<pair<int, double>> s0(1, make_pair(0, one));
  vector<vector<pair<int, double>>> s1(1, s0);
  lmtuv_.push_back(s1);

  for (int l = 1; l != LEND; ++l) {

    map<int, int> mapping;
    int cnt = 0;
    for (int z = 0; z <= l; ++z) {
      for (int y = 0; y <= l - z; ++y) {
        const int x = l - y - z;
        if (x < 0) continue;
        const int key = x + y * LARGE + z * LARGE * LARGE;
        mapping.insert(make_pair(key, cnt));
        ++cnt;
      }
    }

    vector<vector<pair<int, double>>> mtuv;
    for (int n = 0; n != 2 * l + 1; ++n) {
      int m = l - (n / 2);
      if (n % 2 == 1) m *= -1;

      const int vm2 = m < 0 ? 1 : 0;
      const int absm = m > 0 ? m : -m;
      vector<pair<int, double>> tuv;

      const double Nlms = one / pow(two, absm) / factorial[l] * sqrt((m == 0 ? one : two) * factorial[l + absm] * factorial[l - absm]);
      const int tmax = floor((l - absm) / 2.0);
      for (int t = 0; t <= tmax; ++t) {
        for (int u = 0; u <= t; ++u) {
          const int vmax = 2 * floor((absm - vm2) / 2.0) + vm2;
          for (int v2 = vm2; v2 <= vmax; v2 += 2) {
            assert((v2 - vm2) % 2 == 0);
            const double Clmtuv = pow(-one, t + ((v2 - vm2) / 2)) * pow(quarter, t)
                                * data.comb(l, t) * data.comb(l - t, absm + t)
                                * data.comb(t, u) * data.comb(absm, v2) * Nlms;
            const int xexp = 2 * t + absm - 2 * u - v2;
            const int yexp = 2 * u + v2;
            const int zexp = l - 2 * t - absm;
            double denom = one;
            map<int, int>::const_iterator current = mapping.find(xexp + yexp * LARGE + zexp * LARGE * LARGE);
            assert(current != mapping.end());
            const double coeff = Clmtuv / denom;
            tuv.push_back(make_pair(current->second, coeff));
          }
        }
      }
      /* I have no idea if this will work for the general case, but it works up to d */
      double scale = 0.0;
      for(auto ituv = tuv.begin(); ituv != tuv.end(); ++ituv) {
        scale += ituv->second*ituv->second;
      }

      for(auto ituv = tuv.begin(); ituv != tuv.end(); ++ituv) {
        ituv->second /= scale;
      }

      mtuv.push_back(tuv);
    }


    lmtuv_.push_back(mtuv);
  }
}

vector<double> MoldenIn::transform_cart(vector<double> carts, int ang_l) {
   vector<vector<pair<int,double>>> mtuv = lmtuv_.at(ang_l);

   vector<double> out;
   for(auto& im : mtuv) {
     double value = 0.0;

     for(auto& ituv : im) {
        value += (ituv.second) * carts.at(ituv.first);
     }

     out.push_back(value);
   }

   return out;
}
