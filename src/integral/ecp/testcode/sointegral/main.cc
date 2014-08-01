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

  cout << " ++++ TEST int r^n exp(-zeta r*r) <phi_A|lm><lm|l|lm'><lm'|phi_C> ++++ " << endl;
  cout << endl;
  cout << " ++++++++ INITIALISATION ++++++++ " << endl;
  const array<double, 3> centreB = {{0.0000, 0.0000, 0.0000}};
  const array<int, 2> lm = {{1, 0}};
  auto rsh = make_shared<const SphHarmonics>(lm, centreB);
  rsh->print();

  const array<double, 3> centreA = {{0.0000, 0.0000, 0.0000}};
  const array<int, 3> angular_momentumA = {{0, 1, 1}};
  const double alphaA = 1.0;
  auto cargaussA = make_shared<const CartesianGauss>(alphaA, angular_momentumA, centreA);
  cargaussA->print();

  const array<double, 3> centreC = {{0.0000, 0.0000, 0.0000}};
  const array<int, 3> angular_momentumC = {{0 , 1, 1}};
  const double alphaC = 1.0;
  auto cargaussC = make_shared<const CartesianGauss>(alphaC, angular_momentumC, centreC);
  cargaussC->print();

  cout << "Using unnormalized gaussian..." << endl;


  const int max_iter = 100;
  const double thresh_int = 10e-5;
  const double zeta = 0.0;
  const int n = 2;
  const int ic = 1;

  tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>, shared_ptr<const SphHarmonics>, double, int, int>
    so(cargaussA, cargaussC, rsh, zeta, n, ic);

  SOIntegral soint(so);
  const double r = 1.0;
  const double f = soint.compute(r);
  cout << setw(20) << " r = " << setprecision(12) << r
       << setw(20) << " f(r) = " << setprecision(12) << f << endl;


  RadialInt<SOIntegral, tuple<shared_ptr<const CartesianGauss>, shared_ptr<const CartesianGauss>,
                        shared_ptr<const SphHarmonics>, double, int, int>> rad(so, true, max_iter, thresh_int);

#if 0
  /* CHECKED */
  cout << " Test Radial Integration " << endl;
  RadialInt<Gaussian_Int, const double> radial(zeta, true, max_iter, thresh_int);
  cout << "Analytic = " << std::sqrt(pi__ / zeta) / 2.0 << endl;
#endif

  return 0;

}

