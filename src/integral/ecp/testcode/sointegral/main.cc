//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 28, 2014
//

#include <iomanip>
#include "radialquad.h"
#include "so.h"

using namespace std;
using namespace bagel;

class Gaussian_Int {
  protected:
    double exponent_;

  public:
    Gaussian_Int(const double exp) : exponent_(exp) {}
    ~Gaussian_Int() {}

    double compute(const double r) { return exp(-exponent_ * r * r); }
    double exponent() { return exponent_; }
};

int main() {

  /* ++++ TEST int r^n exp(-zeta r*r) <phi_A|lm><lm|l|lm'><lm'|phi_C> ++++  */
  cout << endl;
  const array<double, 3> centreA = {{0.0000, 0.0000, 0.0000}};
  const array<int, 3> angular_momentumA = {{0, 1, 0}};
  const double alphaA = 1.0;
  auto cargaussA = make_shared<const CartesianGauss>(alphaA, angular_momentumA, centreA);
  cargaussA->print();

  const array<double, 3> centreC = {{0.0000, 0.0000, 0.0000}};
  const array<int, 3> angular_momentumC = {{1 , 0, 0}};
  const double alphaC = 1.0;
  auto cargaussC = make_shared<const CartesianGauss>(alphaC, angular_momentumC, centreC);
  cargaussC->print();

  cout << "Using unnormalized gaussian..." << endl;


  const int max_iter = 100;
  const double thresh_int = 10e-5;

  /* ECP Parameters */
  const double ecp_coef = 1.0;
  const double ecp_exp = 0.0;
  const int ecp_r = 0;

  for (int ic = 0; ic != 3; ++ic) {
    double out = 0.0;
    const array<double, 3> centreB = {{0.0000, 0.0000, 0.0000}};
    const int l = 1;
    const array<int, 2> lm = {{l, 0}}; // m can be any number st |m| <= l
    auto rsh = make_shared<const SphHarmonics>(lm, centreB);
    tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>, shared_ptr<const SphHarmonics>, double, int, int>
      so(cargaussA, cargaussC, rsh, ecp_exp, ecp_r, ic);
    //rsh->print();
    //SOIntegral soint(so);
    //soint.compute(1.461297473945);
    RadialInt<SOIntegral, tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>,
              shared_ptr<const SphHarmonics>, double, int, int>> rad(so, false, max_iter, thresh_int);
    out = ecp_coef*rad.integral();
    cout << "int r^n exp(-zrr) <" << angular_momentumA[0] << angular_momentumA[1] << angular_momentumA[2]
         << "|" << l << "m><" << l << "m|l_" << ic << "|" << l << "m'><" << l << "m'|"
         << angular_momentumC[0] << angular_momentumC[1] << angular_momentumC[2] << "> dr" << endl;
    cout << "out[" << ic << "] = " << setprecision(12) << out << endl;
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

