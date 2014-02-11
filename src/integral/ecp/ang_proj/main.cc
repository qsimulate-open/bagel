//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#include "proj.h"

using namespace std;

int main() {


  cout << " Expansion of a Gaussian about a different centre " << endl;
  BesselI besselI;
  const int l = 3;
  const double x =0.1;
  cout << "I_n(x) where n = " << l << " and x = " << x << " is   " << besselI.besseln(l, x) << endl;

  SH sh;
  const int m = 3;
  cout << "(l m x) = (" << l << "  " << m << "  " << x << ")   P_lm(x) =  " << sh.alegendre(l, fabs(m), x) << endl;
  const double phi = 1.3;
  const double theta = 3.1;
  cout << "th = " << theta << "  ph = " << phi << "  Y_lm(x) = " << sh.ylm(l, m, theta, phi).real() << "    " << sh.ylm(l, m, theta, phi).imag() << endl;

  GaussOntoSph gos;
  cout << "coef =   " << gos.compute_c(2, 0, 2, 0, 0) << endl;
  std::list<std::shared_ptr<CartesianGauss>> gauss;

  gauss = gos.sphcar(1, 1);
  cout << "No. of Cartesian Gaussians = " << gauss.size() << endl;
  for (auto& it : gauss) {
    cout << setw(17) << setprecision(9) << it->angular_momentum()[0] << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum()[1] << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum()[2] << endl;
  }

  return 0;

}
