//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 2015
//

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <gsl/gsl_sf_gamma.h>
#include <src/util/math/gamma.h>

using namespace std;

int main() {

// Test Gamma function
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

  return 0;

}

