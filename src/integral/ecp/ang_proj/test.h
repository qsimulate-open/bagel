//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#ifndef __SRC_INTEGRAL_ECP_ANG_PROJ_TEST_H
#define __SRC_INTEGRAL_ECP_ANG_PROJ_TEST_H

#include "proj.h"

using namespace mpfr;
using namespace boost;

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
      return static_cast<mpreal>(int1 * int2);
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

#endif
