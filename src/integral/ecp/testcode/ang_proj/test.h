//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#ifndef __SRC_INTEGRAL_ECP_ANG_PROJ_TEST_H
#define __SRC_INTEGRAL_ECP_ANG_PROJ_TEST_H

#include <cassert>
#include "proj.h"

using namespace mpfr;
using namespace boost;

namespace test {

class Gaussian_Int {

  protected:
    double exponent_;

  public:
    Gaussian_Int(const double exp) : exponent_(exp) {}
    ~Gaussian_Int() {}

    mpreal compute(const mpreal r) { return exp(static_cast<mpreal>(-exponent_) * r * r); }

    double exponent() { return exponent_; }

};

class Projection2 {

  protected:
    std::shared_ptr<ProjectionInt> projAB_;
    std::shared_ptr<ProjectionInt> projCB_;

  public:

    Projection2(std::pair<std::shared_ptr<ProjectionInt>, std::shared_ptr<ProjectionInt>> projs) :
      projAB_(projs.first),
      projCB_(projs.second)
    {}
    ~Projection2() {}

    mpreal compute(const mpreal r) {
      const double int1 = projAB_->compute(r.toDouble());
      const double int2 = projCB_->compute(r.toDouble());
      return static_cast<mpreal>(int1 * int2 * r * r);
    }
};

class GaussianProduct {

  protected:
    int iX_;
    std::shared_ptr<CartesianGauss> gA_;
    std::shared_ptr<CartesianGauss> gB_;

  public:
    GaussianProduct(std::tuple<int, std::shared_ptr<CartesianGauss>, std::shared_ptr<CartesianGauss>> gproduct) :
      iX_(std::get<0>(gproduct)),
      gA_(std::get<1>(gproduct)),
      gB_(std::get<2>(gproduct))
    {}
    ~GaussianProduct() {}

    mpreal compute(const mpreal r) {
      if (r.toDouble() > 15.0) {
        return static_cast<mpreal>(0.0);
      } else {
        const double rd = r.toDouble();
        const double g1 = gA_->compute1D(iX_, rd + gA_->centre(iX_));
        const double g2 = gB_->compute1D(iX_, rd + gB_->centre(iX_));
        return static_cast<mpreal>(g1 * g2);
      }
    }

};

class ABBA_ss {

  protected:
    std::shared_ptr<CartesianGauss> gA_;
    std::shared_ptr<RealSH> shB_;

  public:
     ABBA_ss(std::pair<std::shared_ptr<CartesianGauss>, std::shared_ptr<RealSH>> gsh) :
       gA_(gsh.first),
       shB_(gsh.second)
     {}
     ~ABBA_ss() {}

     mpreal compute(const mpreal r) {
       const double pi = static_cast<double>(atan(1.0) * 4.0);
       std::array<double, 3> AB;
       for (int i = 0; i != 3; ++i) AB[i] = shB_->centre(i) - gA_->centre(i);
       const double dABsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
       const double x = 2.0 * gA_->exponent() * std::sqrt(dABsq) * r.toDouble();
       std::shared_ptr<Modified_Spherical_Bessel_Iexp> msbessel = std::make_shared<Modified_Spherical_Bessel_Iexp>(0);
       const mpreal sbessel = msbessel->compute(x);
//     const double sbessel = boost::math::sph_bessel(0, x);
//     const mpreal sbessel_0 = static_cast<mpreal>((exp(x) - exp(-x)) / (2.0 * x));
//     std::cout << "sbessel = " << sbessel << std::endl;
//     std::cout << "sbessel_0 = " << sbessel_0.toDouble() << std::endl;
//     assert(sbessel == sbessel_0.toDouble());
       const double prefactor = 4.0 * pi * std::pow(2.0 * gA_->exponent() / pi, 1.5);
       return static_cast<mpreal>(prefactor * sbessel * sbessel * r * r * exp(-2.0 * (gA_->exponent() - r) * (gA_->exponent() - r)));
     }

};

class ABBB_ss {

  protected:
    std::shared_ptr<CartesianGauss> gA_;
    std::shared_ptr<RealSH> shB_;

  public:
     ABBB_ss(std::pair<std::shared_ptr<CartesianGauss>, std::shared_ptr<RealSH>> gsh) :
       gA_(gsh.first),
       shB_(gsh.second)
     {}
     ~ABBB_ss() {}

     mpreal compute(const mpreal r) {
       const double pi = static_cast<double>(atan(1.0) * 4.0);
       std::array<double, 3> AB;
       for (int i = 0; i != 3; ++i) AB[i] = shB_->centre(i) - gA_->centre(i);
       const double dABsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
       const double x = 2.0 * gA_->exponent() * std::sqrt(dABsq) * r.toDouble();
       std::shared_ptr<Modified_Spherical_Bessel_Iexp> msbessel = std::make_shared<Modified_Spherical_Bessel_Iexp>(0);
       const mpreal sbessel = msbessel->compute(x);
       const double prefactor = 4.0 * pi * std::pow(2.0 * gA_->exponent() / pi, 1.5);
       return static_cast<mpreal>(prefactor * sbessel * r * r * exp(-r * r) * exp(-(gA_->exponent() - r) * (gA_->exponent() - r)));
     }

};

}

#endif
