//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#include "proj.h"
//#include "src/integral/carsphlist.h"

using namespace std;

int main() {

#if 0
  cout << " Expansion of a Gaussian about a different centre " << endl;
  BesselI besselI;
  int l = 0;
  const double x =0.1;
  cout << "I_n(x) where n = " << l << " and x = " << x << " is   " << besselI.besseln(l, x) << endl;

  SH sh;
  const int m = 0;
  cout << "(l m x) = (" << l << "  " << m << "  " << x << ")   P_lm(x) =  " << sh.alegendre(l, fabs(m), x) << endl;
  const double phi = 1.3;
  const double theta = 3.1;
  cout << "th = " << theta << "  ph = " << phi << "  Y_lm(x) = " << sh.ylm(l, m, theta, phi).real() << "    " << sh.ylm(l, m, theta, phi).imag() << endl;
#endif

#if 0
  GaussOntoSph gos;
  cout << "coef =   " << gos.compute_c(2, 0, 2, 0, 0) << endl;

  std::list<std::shared_ptr<CartesianGauss>> gauss;
  std::array<double, 3> centre = {0.0, 0.0, 0.0};
  gauss = gos.sphcar(centre, 0,0);
  cout << "No. of Cartesian Gaussians = " << gauss.size() << endl;
  for (auto& it : gauss) {
    cout << setw(17) << setprecision(9) << it->angular_momentum(0) << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum(1) << " ";
    cout << setw(17) << setprecision(9) << it->angular_momentum(2) << endl;
  }
  for (auto& it : gauss) {
    cout << "Centre of Gaussian =  " << endl;
    for (int i = 0; i != 3; ++i) {
      cout << setw(17) << setprecision(9) << it->centre(i) << "   ";
    }
    cout << endl;
  }
#endif

#if 0
const int maxl = 3;
for (int iz = 0; iz <= maxl; ++iz) {
  for (int iy = 0; iy <= maxl - iz; ++iy) {
    const int ix = maxl - iz - iy;
    Comb comb;
    std::array<int, 3> cartesian = {ix, iy, iz};
    std::shared_ptr<CarSph> carsph = std::make_shared<CarSph>(cartesian);
    carsph->transform_CarSph();
    carsph->print();
    cout << "---" << endl;
  }
}
#endif

#if 1
  std::array<double, 3> centreB = {0.0, 0.0, 0.0};
  std::array<int, 2> lm = {1, 0};
  std::shared_ptr<RealSH> rsh = std::make_shared<RealSH>(lm, centreB);
  rsh->print();

  std::array<double, 3> centreA = {0.0, 0.0, 0.0};
  std::array<int, 3> angular_momentum = {0, 0, 1};
  const double alpha = 0.0;
  std::shared_ptr<CartesianGauss> cargauss = std::make_shared<CartesianGauss>(alpha, angular_momentum, centreA);
  cargauss->print();

  AngularProj proj(cargauss, rsh);
  cout << "Now integrate... " << endl;
  const double integral = proj.integrate(1.0);
  cout << " Ans = " << integral << endl;

  cout << endl;

  const int lp = 3;
  for (int i = 0; i <= 2*lp; ++i) {
    const int mp = i - lp;
    std::array<int, 2> lmp = {lp, mp};
    std::shared_ptr<SphUSP> sphusp = std::make_shared<SphUSP>(lmp);

//  sphusp->print();

  }
#endif
  cout << "***  TEST INTEGRATION ***" << endl;
  std::array<int, 3> ijk = {1, 4, 2};
  std::pair<int, int> lm1 = std::make_pair(4, 2);
  std::pair<int, int> lm2 = std::make_pair(3, 1);
  const double ans = proj.integrate2SH1USP(lm1, lm2, ijk);
  cout << "int_(lm1 * lm2 * xyz) = " << ans << endl;


  return 0;

}
