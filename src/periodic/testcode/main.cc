//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 2015
//

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <gsl/gsl_sf_gamma.h>
#include <src/util/math/gamma.h>
#include <src/periodic/testcode/compute_mlm.h>

using namespace std;

void test_gamma_upper() {

  // Test upper incomplete gamma function
  double x = 1e-3;
  bagel::Gamma_upper gamma;
  for (int i = 0; i != 100; ++i) {
    for (int j = 0; j != 10; ++j) {
      const double l = j + 0.5;
      cout << "(l, x) = (" << l << ", " << setprecision(5) << x << ")" << endl;
      cout << "boost::tgamma                  = " << setw(20) << setprecision(12) << boost::math::tgamma(l, x) << "  "
          << "erfc(x)                        = " << setw(20) << setprecision(12) << boost::math::erfc(x) << endl;
      cout << "bagel::Gamma_upper             = " << setw(20) << setprecision(12) << gamma(j, x) << "  "
          << "erfc(x)                        = " << setw(20) << setprecision(12) << erfc(x) << endl;
      cout << "gsl_sf_gamma_inc               = " << setw(20) << setprecision(12) << gsl_sf_gamma_inc(l, x) << endl;
      cout << endl;
    }
    x += 0.01;
  }

};


void test_gamma_lower_scaled() {

  // Test scaled lower gamma function for large argument
  double x = 100.0;
  const double beta = sqrt(pi__);
  bagel::Gamma_lower_scaled gamma_s;
  for (int i = 0; i != 10; ++i) {
    for (int j = 0; j != 10; ++j) {
      const double l = j + 0.5;
      const double z = pi__ * pow(x, 2.0);
      const double boost_gamma = boost::math::tgamma_lower(l, z)/pow(x, j+1.0);
      const double bagel_gamma = gamma_s(j, z, beta);
      const double gsl_gamma = (gsl_sf_gamma(l) - gsl_sf_gamma_inc(l, z))/pow(x, j+1.0);
      if (abs(bagel_gamma - boost_gamma) > 1e-14 || abs(bagel_gamma - gsl_gamma) > 1e-14 || abs(boost_gamma - gsl_gamma) > 1e-14) {
        cout << "(l, z) = (" << l << ", " << setprecision(5) << z << ")" << endl;
        cout << "boost      = " << setw(20) << setprecision(12) << boost_gamma << endl;
        cout << "bagel      = " << setw(20) << setprecision(12) << bagel_gamma << endl;
        cout << "gsl        = " << setw(20) << setprecision(12) << gsl_gamma << endl << endl;
      }
    }
    x += 0.1;
  }

};
    if (abs(bagel_gamma - boost_gamma) > 1e-14 || abs(bagel_gamma - gsl_gamma) > 1e-14 || abs(boost_gamma - gsl_gamma) > 1e-14) {
      cout << "(l, z) = (" << l << ", " << setprecision(5) << z << ")" << endl;
      cout << "boost      = " << setw(20) << setprecision(12) << boost_gamma << endl;
      cout << "bagel      = " << setw(20) << setprecision(12) << bagel_gamma << endl;
      cout << "gsl        = " << setw(20) << setprecision(12) << gsl_gamma << endl << endl;
    }
  }
  x += 0.1;
}

};


void test_mlm() {

  vector<complex<double>> mlm;
  const int ws = 1;
  const int lmax = 10;
  const int limit = 10;
  const int thresh = PRIM_SCREEN_THRESH;
  mlm = compute_mlm(ws, lmax, limit, thresh);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int mm = 0; mm <= 2*l+1; ++mm, ++cnt) {
      const int m = mm - l;
      if (abs(mlm[cnt].real()) > 1e-10)
        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << mlm[cnt].real() << endl;
    }
  }

};


int main() {

//  test_gamma_upper();
//  test_mlm();
  test_gamma_lower_scaled();
  return 0;
}
