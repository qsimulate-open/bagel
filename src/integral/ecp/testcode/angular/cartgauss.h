//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 28, 2014
//


#ifndef __ECP_TESTCODE_SOINTEGRAL_CARTGAUSS_H
#define __ECP_TESTCODE_SOINTEGRAL_CARTGAUSS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <memory>
#include "src/util/constants.h"
#include "src/util/math/factorial.h"

using namespace bagel;
using namespace std;

const static DoubleFactorial dfact;
const static Factorial fact;

namespace test {

class CartesianGauss {
  protected:

    double norm_;
    double exponent_;
    double weight_;
    array<int, 3> angular_momentum_;
    array<double, 3> centre_;

  public:

    CartesianGauss(const double alpha, const array<int, 3> l, const array<double, 3> c) :
      exponent_(alpha),
      angular_momentum_(l),
      centre_(c)
      { norm_ = normalise(); }
    ~CartesianGauss() {}

    double normalise() const {
      Factorial fact;
      const int l = angular_momentum_[0] + angular_momentum_[1] + angular_momentum_[2];
      const double fact1 = fact(angular_momentum_[0]) * fact(angular_momentum_[1]) * fact(angular_momentum_[2]);
      const double fact2 = fact(2*angular_momentum_[0]) * fact(2*angular_momentum_[1]) * fact(2*angular_momentum_[2]);

      return pow(2 * exponent_ / pi__, 0.75) * sqrt(pow(4 * exponent_, l) * pow(2, l) * fact1 / fact2);
    }

    int angular_momentum(const int i) const { return angular_momentum_[i]; }
    array<int, 3> angular_momentum() const { return angular_momentum_; }
    int total_ang() const { return angular_momentum_[0]+angular_momentum_[1]+angular_momentum_[2]; }

    double centre(const int i) const { return centre_[i]; }
    array<double, 3> centre() const { return centre_; }

    double exponent () const { return exponent_; }

    void set_weight(const double wt) { weight_ = wt; }

    double compute(const array<double, 3> centre) {
      array<double, 3> RA;
      for (int i = 0; i != 3; ++i) RA[i] = centre[i] - centre_[i];
      double rsq = RA[0] * RA[0] + RA[1] * RA[1] + RA[2] * RA[2];
      return norm_ * pow(RA[0], angular_momentum_[0]) * pow(RA[1], angular_momentum_[1]) * pow(RA[2], angular_momentum_[2]) * exp(-exponent_ * rsq);
    }

    double compute1D(const int iX, const double centreX) {
      Factorial fact;
      const int lx = angular_momentum_[iX];
      const double normX = pow(2 * exponent_ / pi__, 0.25) * sqrt(pow(4 * exponent_, lx) * pow(2, lx) * fact(lx) / fact(2*lx));
      double r = centreX - centre_[iX];
      return normX * pow(r, lx) * exp(-exponent_ * r * r);
    }

    void print() const {
      cout << "** Cartesian Gaussian **" << endl;
      cout << "Angular momentum (lx, ly, lz) = (" << angular_momentum_[0] << ", "
                                                       << angular_momentum_[1] << ", "
                                                       << angular_momentum_[2] << ")" << endl;
      cout << "Centre = " ;
      for (int i = 0; i != 3; ++i) {
        cout << setw(17) << setprecision(9) << centre_[i] << ";   ";
      }
      cout << endl;
      cout << "Exponent alpha = " << setw(17) << setprecision(9) << exponent_  << endl;
      cout << "Normalisation const N = " << setw(17) << setprecision(9) << norm_ << endl;
      cout << "---" << endl;
    }

};

}

#endif
