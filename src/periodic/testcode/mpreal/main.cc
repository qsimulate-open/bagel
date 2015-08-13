//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 2015
//

#include <boost/math/special_functions/erf.hpp>
#include <gsl/gsl_sf_gamma.h>
//#include <src/util/math/gamma.h>
#include <src/periodic/testcode/compute_mlm.h>

using namespace std;

mpreal half = "0.5";
mpreal two  = "2.0";

void test_gamma_upper() {

  // Test upper incomplete gamma function
  double dx = 1e-13;
  for (int i = 0; i != 100; ++i) {
    mpreal x = static_cast<mpreal>(dx);
    for (int j = 0; j != 10; ++j) {
      const mpreal l = j + half;
      cout << "(l, x) = (" << l << ", " << setprecision(5) << x << ")" << endl;
      cout << "boost::tgamma                  = " << setw(20) << setprecision(12) << boost::math::tgamma(l.toDouble(), dx) << endl;
      mpreal bagel_gamma = compute_gamma_upper(j, x);
      cout << "bagel::Gamma_upper             = " << setw(20) << setprecision(12) << bagel_gamma.toDouble() << endl;
      cout << "gsl_sf_gamma_inc               = " << setw(20) << setprecision(12) << gsl_sf_gamma_inc(l.toDouble(), dx) << endl;
      cout << endl;
    }
    dx += 0.01;
  }

};


void test_gamma_lower_scaled() {

  // Test scaled lower gamma function for large argument
  double dx = 100.0;
  mpreal beta = GMPPISQRT;
  mpreal pi = GMPPI;
  for (int i = 0; i != 10; ++i) {
    mpreal x = static_cast<mpreal>(dx);
    for (int j = 0; j != 10; ++j) {
      const mpreal l = j + half;
      mpreal z = beta * beta * pow(x, two);
      const double boost_gamma = boost::math::tgamma_lower(l.toDouble(), z.toDouble())/std::pow(dx, j+1.0);
      mpreal bagel_gamma = compute_gamma_lower_scaled(j, z, beta);
      const double gsl_gamma = (gsl_sf_gamma(l.toDouble()) - gsl_sf_gamma_inc(l.toDouble(), z.toDouble()))/std::pow(dx, j+1.0);
      if (abs(bagel_gamma.toDouble() - boost_gamma) > 1e-14 || abs(bagel_gamma.toDouble() - gsl_gamma) > 1e-14 || abs(boost_gamma - gsl_gamma) > 1e-14) {
        cout << "(l, z) = (" << j << ", " << setprecision(5) << z << ")" << endl;
        cout << "boost      = " << setw(20) << setprecision(12) << boost_gamma << endl;
        cout << "bagel      = " << setw(20) << setprecision(12) << bagel_gamma.toDouble() << endl;
        cout << "gsl        = " << setw(20) << setprecision(12) << gsl_gamma << endl << endl;
      }
    }
    dx += 0.1;
  }

};


void test_gamma() { //G(l+1/2)

  // Test complete gamma function for half integers
  for (int l = 150; l != 300; ++l) {
    const int n = 2 * l + 1;
    const double boost_gamma = boost::math::tgamma(n/2.0);
    const mpreal bagel_gamma = compute_gamma(n);
    const double gsl_gamma = gsl_sf_gamma(n/2.0);
    if (abs(bagel_gamma.toDouble() - boost_gamma) > 1e-14 || abs(bagel_gamma.toDouble() - gsl_gamma) > 1e-14 || abs(boost_gamma - gsl_gamma) > 1e-14) {
      cout << "l = " << l << endl;
      cout << "boost      = " << setw(20) << setprecision(12) << boost_gamma << endl;
      cout << "bagel      = " << setw(20) << setprecision(12) << bagel_gamma.toDouble() << endl;
      cout << "gsl        = " << setw(20) << setprecision(12) << gsl_gamma << endl << endl;
    }
  }

};


void test_mlm() {

  vector<complex<mpreal>> mlm;
  const int ws = 1;
  const int lmax = 5;
  const int limit = 10;
  const mpreal thresh = static_cast<mpreal>(PRIM_SCREEN_THRESH);
  mlm = compute_mlm(ws, lmax, limit, thresh);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int mm = 0; mm <= 2*l+1; ++mm, ++cnt) {
      const int m = mm - l;
      if (abs((mlm[cnt].real()).toDouble()) > 1e-10)
        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << (mlm[cnt].real()).toDouble() << endl;
    }
  }

};


int main() {

//  test_gamma_upper();
  test_mlm();
//  test_gamma_lower_scaled();
//  test_gamma();

  return 0;
}
