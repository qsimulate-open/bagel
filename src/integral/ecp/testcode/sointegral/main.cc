//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 28, 2014
//

#include "radialquad.h"
#include "so.h"

using namespace std;
using namespace bagel;

int main() {

  cout << " ++++ TEST <phi_A|lm><lm|l|lm'><lm'|phi_C> ++++ " << endl;
  cout << endl;
  cout << " ++++++++ INITIALISATION ++++++++ " << endl;
  const array<double, 3> centreB = {{0.0000, 0.0000, 0.0000}};
  const array<int, 2> lm0 = {{0, 0}};
  auto rsh0 = make_shared<SphHarmonics>(lm0, centreB);
  rsh0->print();

  const array<int, 2> lm1 = {{0, 0}};
  auto rsh1 = make_shared<SphHarmonics>(lm1, centreB);
  rsh1->print();

  const array<double, 3> centreA = {{0.0000, 0.0000, 0.5000}};
  const array<int, 3> angular_momentumA = {{0, 2, 0}};
  const double alphaA = 1.0;
  auto cargaussA = make_shared<CartesianGauss>(alphaA, angular_momentumA, centreA);
  cargaussA->print();

  const array<double, 3> centreC = {{0.0000, 0.0000, 0.0000}};
  const array<int, 3> angular_momentumC = {{2, 0, 0}};
  const double alphaC = 1.0;
  auto cargaussC = make_shared<CartesianGauss>(alphaC, angular_momentumC, centreC);
  cargaussC->print();

  cout << "Using unnormalized gaussian..." << endl;

  return 0;

}

