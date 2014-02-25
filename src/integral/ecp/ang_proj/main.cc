//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#include "proj.h"
#include "src/integral/carsphlist.h"

using namespace std;

int main() {


  cout << " Expansion of a Gaussian about a different centre " << endl;
  BesselI besselI;
  int l = 1;
  const double x =0.1;
  cout << "I_n(x) where n = " << l << " and x = " << x << " is   " << besselI.besseln(l, x) << endl;

  SH sh;
  const int m = 0;
  cout << "(l m x) = (" << l << "  " << m << "  " << x << ")   P_lm(x) =  " << sh.alegendre(l, fabs(m), x) << endl;
  const double phi = 1.3;
  const double theta = 3.1;
  cout << "th = " << theta << "  ph = " << phi << "  Y_lm(x) = " << sh.ylm(l, m, theta, phi).real() << "    " << sh.ylm(l, m, theta, phi).imag() << endl;

  GaussOntoSph gos;
  cout << "coef =   " << gos.compute_c(2, 0, 2, 0, 0) << endl;
  std::list<std::shared_ptr<CartesianGauss>> gauss;

  std::array<double, 3> centre = {1.0, 1.0, 1.0};
  gauss = gos.sphcar(centre, 1,-1);
  cout << "No. of Cartesian Gaussians = " << gauss.size() << endl;
  for (auto& it : gauss) {
    cout << setw(17) << setprecision(9) << it->angular_momentum(0) << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum(1) << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum(2) << endl;
  }
  for (auto& it : gauss) {
    for (int i = 0; i != 3; ++i) {
      cout << "Centre = " << it->centre(i) << endl;
    }
  }

  Comb comb;
  std::array<int, 3> cartesian = {0, 0, 2};
  std::shared_ptr<CarSph> carsph = std::make_shared<CarSph>(cartesian);
  carsph->transform_CarSph();
  carsph->print();

  return 0;

}
