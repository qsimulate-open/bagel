//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#include "radial.h"
#include "test.h"
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
  std::array<double, 3> centreB = {0.0, 0.0, 1.420616};
  std::array<int, 2> lm = {0, 0};
  std::shared_ptr<RealSH> rsh = std::make_shared<RealSH>(lm, centreB);
  rsh->print();

  std::array<double, 3> centreA = {0.0, 0.0, 0.305956};
  std::array<int, 3> angular_momentumA = {0, 0, 0};
  const double alphaA = 1.0;
  std::shared_ptr<CartesianGauss> cargaussA = std::make_shared<CartesianGauss>(alphaA, angular_momentumA, centreA);
  cargaussA->print();

  std::array<double, 3> centreC = {0.0, 0.0, 0.305956};
  std::array<int, 3> angular_momentumC = {0, 0, 0};
  const double alphaC = 1.0;
  std::shared_ptr<CartesianGauss> cargaussC = std::make_shared<CartesianGauss>(alphaC, angular_momentumC, centreC);
  cargaussC->print();

  const double r = 1.0;
  cout << " r  = " << r << endl;

  std::shared_ptr<ProjectionInt> projAB = std::make_shared<ProjectionInt>(cargaussA, rsh);
  const double int1 = projAB->compute(r);
  cout << " < phi_A | lm_B >(r)  =  " << int1 << endl;

  std::shared_ptr<ProjectionInt> projCB = std::make_shared<ProjectionInt>(cargaussA, rsh);
  const double int2 = projCB->compute(r);
  cout << " < phi_C | lm_B >(r)  =  " << int2 << endl;

  cout << endl;
  cout << " < phi_A | lm_B > < lm_B | phi_C >(r) =  " << int1 * int2 << endl;
#endif

#if 0
  const int lp = 3;
  for (int i = 0; i <= 2*lp; ++i) {
    const int mp = i - lp;
    std::array<int, 2> lmp = {lp, mp};
    std::shared_ptr<SphUSP> sphusp = std::make_shared<SphUSP>(lmp);

//  sphusp->print();

  }
  cout << "***  TEST INTEGRATION ***" << endl;
  std::array<int, 3> ijk = {1, 4, 2};
  std::pair<int, int> lm1 = std::make_pair(4, 2);
  std::pair<int, int> lm2 = std::make_pair(3, 1);
  const double ans = projAB.integrate2SH1USP(lm1, lm2, ijk);
  cout << "int_(lm1 * lm2 * xyz) = " << ans << endl;
#endif

  cout << " Test Radial Integration " << endl;
  const int max_iter = 100;
  const double thresh_int = 10e-10;
  const double exp = 2.0;
  const double pi = static_cast<double>(atan(1.0) * 4.0);

  Radial_Int<Gaussian_Int, const double> radial(max_iter, thresh_int, exp);
  cout << "Analytic = " << std::sqrt(pi / exp) / 2.0 << endl;

#if 0
  cout << " --- " << endl;
  std::pair<std::shared_ptr<ProjectionInt>, std::shared_ptr<ProjectionInt>> projs(projAB, projCB);
  Radial_Int<Projection2, std::pair<std::shared_ptr<ProjectionInt>, std::shared_ptr<ProjectionInt>>> ecp(max_iter, thresh_int, projs);
#endif

  // check normalization ss
  cout << endl;
  std::pair<std::shared_ptr<CartesianGauss>, std::shared_ptr<CartesianGauss>> gproduct(cargaussA, cargaussA);
  Radial_Int<GaussianProduct, std::pair<std::shared_ptr<CartesianGauss>, std::shared_ptr<CartesianGauss>>> norm(max_iter, thresh_int, gproduct);
  cout << "Should be " << std::pow(2.0 / pi, 1.5) * std::pow(std::sqrt(pi / 2.0) / 2.0, 3.0) << endl;

  return 0;

}

