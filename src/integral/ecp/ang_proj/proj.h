//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#ifndef __SRC_INTEGRAL_ECP_ANG_PROJ_PROJ_H
#define __SRC_INTEGRAL_ECP_ANG_PROJ_PROJ_H

#include <iostream>
#include <iomanip>
#include <list>
#include <cmath>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include "mpreal.h"
#include "src/math/comb.h"
#include <complex>
#include "wigner3j.h"

using namespace mpfr;
using namespace boost;

class Factorial {
  protected:

  public:

    Factorial() {}
    ~Factorial() {}
    const mpreal compute(const int n) {
      mpreal fact = "1.0";
      for (int i = 0; i != n; ++i) {
        const mpreal imp = i;
        fact *= (imp + 1);
      }
      return fact;
    }

};

class CartesianGauss {
  protected:

    double norm_;
    double exponent_;
    double weight_;
    std::array<int, 3> angular_momentum_;

  public:
    explicit CartesianGauss(const double alpha, const std::array<int, 3> l) : exponent_(alpha), angular_momentum_(l) { norm_ = normalise(); }
    ~CartesianGauss() {}

    double normalise() {
      Factorial fact;
      mpreal pi = static_cast<mpreal>(atan(1) * 4);
      mpreal mexp = static_cast<mpreal>(exponent_);
      mpreal mnorm = fact.compute(2*angular_momentum_[0]) * fact.compute(2*angular_momentum_[1]) * fact.compute(2*angular_momentum_[2]) * pow(pi, 1.5);
      int l = angular_momentum_[0] + angular_momentum_[1] + angular_momentum_[2];
      mnorm /= (pow(2, 2*l) * fact.compute(angular_momentum_[0]) * fact.compute(angular_momentum_[1]) * fact.compute(angular_momentum_[2]) * pow(mexp, l+1.5));
      mnorm = static_cast<mpreal>(1.0 / sqrt(mnorm));
      return mnorm.toDouble();
    }

    std::array<int, 3> angular_momentum() { return angular_momentum_; }

    void set_weight(const double wt) { weight_ = wt; }

    double compute(const std::array<double, 3> centre) {
      double rsq = centre[0]*centre[0] + centre[1]*centre[1] * centre[2]*centre[2];
      return norm_ * std::pow(centre[0], angular_momentum_[0]) * std::pow(centre[1], angular_momentum_[1]) * std::pow(centre[2], angular_momentum_[2]) * std::exp(-exponent_*rsq);
    }

};

class GaussOntoSph {
  protected:

  public:
    GaussOntoSph() {}
    ~GaussOntoSph() {}

    double compute_c(const int l, const int m, const int lx, const int ly, const int lz) {
      const int am = fabs(m);
      const int j = static_cast<int>(0.5*(lx + ly - am));
      if (lx + ly - am - 2*j != 0) {
        return 0.0;
      } else {
        Factorial fact;
        mpreal t1, t2, t3;
        const int lmam = l - am;
        t1 = fact.compute(2*lx) * fact.compute(2*ly) * fact.compute(2*lz) * fact.compute(l) * fact.compute(lmam);
        t1 /= (fact.compute(2*l) * fact.compute(lx) * fact.compute(ly) * fact.compute(lz) * fact.compute(l + am));
        t1 = sqrt(t1) / pow(2, l) / fact.compute(l);
        std::cout << "t1 = " << t1 << std::endl;
        Comb comb;
        mpreal zero = "0.0";
        mpreal c = zero;
        for (int i = 0; i != lmam/2 + 1; ++i) {
          if (j < 0 || j > i || i > l) {
            t2 = zero;
          } else {
            t2 = static_cast<mpreal>(comb.c(l, i) * comb.c(i, j) * pow(-1, i) * fact.compute(2*l - 2*i)/ fact.compute(lmam - 2*i));
          }
          std::cout << "t2 = " << t2 << std::endl;
          for (int k = 0; k != j+1; ++k) {
            if (k < j || lx - 2*k < 0 || lx - 2*k < am) {
              t3 = zero;
            } else {
              t3 = static_cast<mpreal>(comb.c(j, k) * comb.c(am, lx - 2*k) * pow(-1, 0.5*(am - lx + 2*k)));
            }
            std::cout << "t3 = " << t3 << std::endl;
            c += t1 * t2 * t3;
          }
        }
        return c.toDouble();
      }
    }

    std::list<std::shared_ptr<CartesianGauss>> sphcar(const int l, const int m) {
      std::list<std::shared_ptr<CartesianGauss>> gauss;
      for (int i = 0; i != l+1; ++i) {
        for (int j = 0; j != l+1; ++j) {
          const int k = l - i - j;
          if (k >= 0) {
            const double c = compute_c(l, m, i, j, k);
            std::cout << "(lm)(ijk) = (" << l << m << ") (" << i << j << k << ")    c = " << c << std::endl;
            if (c > 1e-13) {
              std::array<int, 3> ang = {i, j, k};
              gauss.push_back(std::make_shared<CartesianGauss>(0.0, ang));
            }
          }
        }
      }
      return gauss;
    }

    void angular_integral(const int a, const int b, const int c, const int ld, const int l, const int m) {
       // int{ylm * ypq * x^a * y^b * z^c} = C * int{x^(r+u+a) * y^(s+v+b) * z^(t+w+c)}
       // int{x^i * y^j * z^k} = 0 if i, j, or k odd
       //                      = (i-1)!!(j-1)!!(k-1)!!(i+j+k+1)!! if i, j, or k even
    }

};

class SH {
  protected:
  const mpreal pi_ = "3.1415926535897932384626433";

  public:

  SH() {}
  ~SH() {}

  const double alegendre(const int l, const int m, const double x) {
    if (m < 0 || m > l || fabs(x) > 1.0) throw std::runtime_error("SH: m must be in [0, l] and x in [-1, 1]");
    double pmm = 1.0;
    if (m > 0) {
      double somx2 = std::sqrt((1.0 - x)*(1.0 + x));
      double fact = 1.0;
      for (int i = 0; i != m; ++i) {
        pmm *= -fact * somx2;
        fact += 2.0;
      }
    }
    if (l == m) {
      return pmm;
    } else {
      double pmmp1 = x * (2.0 * m + 1) * pmm;
      if (l == (m+1)) {
        return pmmp1;
      } else {
        double plm = 0.0;
        for (int i = m + 2; i != l; ++i) {
          plm = (x * (2 * i -1) * pmmp1 - (i + m - 1) * pmm) / (i - m);
          pmm = pmmp1;
          pmmp1 = plm;
        }
        return plm;
      }
    }
  }

  const std::complex<double> ylm(const int l, const int m, const double theta, const double phi) {
    Factorial fact;
    const double cth = cos(theta);
    if (fabs(m) > l) throw std::runtime_error ("SH.ylm: |m| > l");

    const int am = fabs(m);
    const double plm = alegendre(l, am, cth);
    mpreal f = fact.compute(l-am)/fact.compute(l+am);
    const double coef = std::sqrt((2*l+1)*f.toDouble()*0.25/pi_.toDouble());
    double real = coef*plm*cos(am*phi);
    double imag = coef*plm*sin(am*phi);
    if (m < 0) {
      real *= std::pow(-1, m);
      imag *= std::pow(-1, m+1);
    }
    std::complex<double> y (real, imag);
    return y;
  }

  const double realSH(const int l, const int m, const double theta, const double phi) {
    Factorial fact;
    const double cth = cos(theta);
    if (fabs(m) > l) throw std::runtime_error ("SH.ylm: |m| > l");

    const int am = fabs(m);
    const double plm = alegendre(l, am, cth);
    mpreal f = fact.compute(l-am)/fact.compute(l+am);
    const double coef = std::sqrt((2*l+1)*f.toDouble()*0.25/pi_.toDouble());
    if (m == 0) {
      return coef*plm;
    } else if (m > 0) {
      return std::pow(-1, m)*std::sqrt(2.0)*coef*plm*cos(am*phi);
    } else {
      return std::pow(-1, m)*std::sqrt(2.0)*coef*plm*sin(am*phi);
    }
  }

};

class RealSH {
  protected:

    std::array<double, 3> centre_;
    int l_;
    int m_;

  public:

};

class BesselI {
  protected:

  public:

    BesselI() {}
    ~BesselI() {}

    const double bessel0(const double x) {
      double bessi;
      double y;
      const double ax = fabs(x);

      if (ax < 3.75) {
        y = x / 3.75;
        y *= y;
        bessi = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
      } else {
        y = 3.75 / ax;
        bessi = (exp(ax) / std::sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 +
               y * (0.916281e-2 + y *(-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
      }

      return bessi;
    }

    const double bessel1(const double x) {
      double bessi;
      double y;
      const double ax = fabs(x);

      if (ax < 3.75) {
        y = x / 3.75;
        y *= y;
        bessi = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
      } else {
        y = 3.75 / ax;
        bessi = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        bessi = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y *(0.163801e-2 + y * (-0.1031555e-1 + y * bessi))));
        bessi *= (exp(ax) / std::sqrt(ax));
      }
      return x < 0.0 ? -bessi : bessi;
    }

    const double besseln(const int l, const double x) {
      const mpreal pi = "3.1415926535897932384626433";
      mpreal mbessi = "0.0";
      Factorial fact;
      const mpreal mfact = fact.compute(l);
      const mpreal half = "0.5";
      if (l == 0) {
        return bessel0(x);
      } else if (l == 1) {
        return bessel1(x);
      } else {
        if (fabs(x) < l * 1e-3) {
          mbessi = static_cast<mpreal>(pow(x * half, l)) / mfact;
          return mbessi.toDouble();
        } else if (x > l * 1e2) {
          mbessi = static_cast<mpreal>(exp(x) / sqrt(2.0 * x * pi));
          return mbessi.toDouble();
        } else {
          if (x == 0) {
            return 0.0;
          } else {
            mpreal tox = static_cast<mpreal>(2.0 / fabs(x));
            mpreal bip = "0.0";
            mpreal bi = "1.0";
            mpreal bim;
            const int k = 2 * (l + std::sqrt(40*l));
            for (int i = k; i > 0; --i) {
              bim = bip + i * tox * bi;
              bip = bi;
              bi = bim;
              if (fabs(bi) > 1e10) {
                mbessi *= 1e-10;
                bi *= 1e-10;
                bip *= 1e-10;
              }
              if (i == l) mbessi = bip;
            }
            mbessi *= static_cast<mpreal>(bessel0(x))/bi;
            const double bessi = mbessi.toDouble();
            return x < 0.0 && (l & 1) ? -bessi : bessi;
          }
        }
      }
    }

};

class grid {

};

class Proj {
  protected:

    std::shared_ptr<const CartesianGauss> gauss_;
    std::shared_ptr<const RealSH> sh_;

  public:

    Proj(std::shared_ptr<const CartesianGauss> gauss, std::shared_ptr<const RealSH> sh) : gauss_(gauss), sh_(sh) {}
    ~Proj() {}
    void compute() {}

};

#endif
