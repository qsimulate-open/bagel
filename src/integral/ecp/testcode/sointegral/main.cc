//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 28, 2014
//

#include <iomanip>
#include "radialquad.h"
#include "so.h"

using namespace std;
using namespace bagel;
using namespace test;

namespace test {

class Gaussian_Int {
  protected:
    double exponent_;

  public:
    Gaussian_Int(const double exp) : exponent_(exp) {}
    ~Gaussian_Int() {}

    double compute(const double r) { return exp(-exponent_ * r * r); }
    double exponent() { return exponent_; }
};

}

using namespace test;
int main() {

  /* ++++ TEST int r^n exp(-zeta r*r) <phi_A|lm><lm|l|lm'><lm'|phi_C> ++++  */
  const array<double, 3> centreA = {{0.0000, 0.0000, 0.0000}};
  const double alphaA = 1.0;
  const int lA = 1;
  const array<double, 3> centreC = {{0.0000, 0.0000, 0.0000}};
  const double alphaC = 1.0;
  const int lC = 1;

  for (int lza = 0; lza <= lA; ++lza)
  for (int lya = 0; lya <= lA-lza; ++lya)
  for (int lzc = 0; lzc <= lC; ++lzc)
  for (int lyc = 0; lyc <= lC-lzc; ++lyc) {
    const int lxa = lA-lza-lya;
    const int lxc = lC-lzc-lyc;
    const array<int, 3> angA = {{lxa, lya, lza}};
    auto cargaussA = make_shared<const CartesianGauss>(alphaA, angA, centreA);
//  cargaussA->print();

    const array<int, 3> angC = {{lxc , lyc, lzc}};
    auto cargaussC = make_shared<const CartesianGauss>(alphaC, angC, centreC);
//  cargaussC->print();

    const array<double, 3> centreB = {{0.0000, 0.0000, 0.0000}};
    const int l = 1;
    const array<int, 2> lm = {{l, 0}}; // m can be any number st |m| <= l
    auto rsh = make_shared<const SphHarmonics>(lm, centreB);
  //rsh->print();
  //SOIntegral soint(so);
  //soint.compute(1.461297473945);

    /* ECP Parameters */
    const double ecp_coef = 1.0;
    const double ecp_exp = 0.0;
    const int ecp_r = 0;

    const int max_iter = 100;
    const double thresh_int = 10e-5;
    /* Using normalized gaussian... */
    cout << "lA = (" << lxa << ", " << lya << ", " << lza << ")    l = " << l << "    ";
    cout << "lC = (" << lxc << ", " << lyc << ", " << lzc << ")    ";
    cout << "(iaa, rab, iab) = (";
    for (int ic = 0; ic != 3; ++ic) {
      tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>, shared_ptr<const SphHarmonics>, double, int, int>
        so(cargaussA, cargaussC, rsh, ecp_exp, ecp_r, ic);
      RadialInt<SOIntegral, tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>,
                shared_ptr<const SphHarmonics>, double, int, int>> rad(so, false, max_iter, thresh_int);
      const double out = ecp_coef*rad.integral();
      if (ic != 2) {
        cout << setprecision(12) << out << ", ";
      } else {
        cout << setprecision(12) << out << ")" << endl;
      }
    }
  }

#if 0
  /* CHECKED */
  cout << " Test Radial Integration " << endl;
  const double zeta = 1.0;
  RadialInt<Gaussian_Int, const double> radial(zeta, true, max_iter, thresh_int);
  cout << "Analytic = " << std::sqrt(pi__ / zeta) / 2.0 << endl;
#endif

  return 0;

}

