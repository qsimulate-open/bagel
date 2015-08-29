//  Author: Hai-Anh Le <anh@u.northwestern.edu>
//  Date: Jul 30, 2014

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <complex>
#include <cassert>

using namespace std;

namespace test {

class Fmlm {
  protected:

  public:
    Fmlm() {}
    ~Fmlm() {}

    double delta (const int i, const int j) const { const double out = (i==j) ? 1.0 : 0.0; return out; }

    complex<double> theta(const int m) const {
      const double tau = (m < 0) ? 0.0 : 1.0;

      const double re = (m==0) ? 0.5 : tau/sqrt(2.0);
      const double im = (m==0) ? 0.0 : (1.0-tau)/sqrt(2.0);

      return complex<double>(re, im);
    }

    array<double, 3> fm0lm1(const int l, const int m0, const int m1) const { /* Im{aa}, Re{ab}, Im{ab} */
      assert(l > 0 && m1 < m0);
      assert(abs(m0) <= l && abs(m1) <=l);
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

};

}

using namespace test;
int main() {

  cout << " +++ TEST <lm|l|lm'> +++ " << endl;

  auto fmm = make_shared<Fmlm>();

  const int lmax = 1;
  for (int l = 0; l <= lmax; ++l) {
    for (int m1 = 0; m1 <= 2*l; ++m1) {
      for (int m2 = 0; m2 <= m1-1; ++m2) {
        const array<double, 3> f = fmm->fm0lm1(l, m1-l, m2-l);
        if (f[0] != 0.0)
          cout << "<" << l << m1-l << "|l_z|" << l << m2-l << "> = " << setprecision(12) << f[0] << endl;
        if (f[1] != 0.0)
          cout << "<" << l << m1-l << "|l_y|" << l << m2-l << "> = " << setprecision(12) << f[1] << endl;
        if (f[2] != 0.0)
          cout << "<" << l << m1-l << "|l_x|" << l << m2-l << "> = " << setprecision(12) << f[2] << endl;
      }
    }
    cout << endl;
  }

  return 0;
}
