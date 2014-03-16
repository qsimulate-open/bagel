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

#endif
