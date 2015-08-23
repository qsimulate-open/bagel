//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#ifndef __SRC_INTEGRAL_ECP_ANG_PROJ_PROJ_H
#define __SRC_INTEGRAL_ECP_ANG_PROJ_PROJ_H

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include "mpreal.h"
#include <src/util/math/comb.h>
#include <src/integral/ecp/wigner3j.h>
#include <src/integral/carsphlist.h>

using namespace mpfr;
using namespace boost;
using namespace bagel;
const static Comb comb;
constexpr int GMPPREC = 256;

namespace test {

class Factorial {
  protected:

  public:

    Factorial() {}
    ~Factorial() {}
    mpreal operator()(const int n) const {
      assert (n >= 0);
      mpreal fact = "1.0";
      for (int i = 0; i != n; ++i) {
        const mpreal imp = i;
        fact *= (imp + 1);
      }
      return fact;
    }

};

class Modified_Spherical_Bessel_Iexp {

  protected:
    int l_;

     mpreal R_l(const double x) const {
       Factorial f;
       mpreal sum = 0.0;
       for (int i = 0; i <= l_; ++i) {
         sum += f(l_ + i) / f(i) / f(l_ -  i) / pow(2.0 * x, i);
       }
       return sum;
     }

  public:

     Modified_Spherical_Bessel_Iexp(const int l) : l_(l) {}
     ~Modified_Spherical_Bessel_Iexp() {}

     mpreal compute(const double x) const {
       const mpreal half = "0.5";
       const mpreal one = "1.0";
       Factorial f;
       if (x < 1e-7) {
         return (one - x) * pow(2.0 * x, l_) * f(l_) / f(2 * l_);
       } else if (x > 16) {
         return half * R_l(-x) / x;
       } else {
         return half * (R_l(-x) - pow(-1.0, l_) * R_l(x) * exp(-2.0 * x)) / x;
       }
     }

};

class SphUSP {
  protected:
    std::array<int, 2> angular_momentum_;
    std::vector<std::pair<double, int>> usp_;

  public:
    SphUSP(const std::array<int, 2> lm) : angular_momentum_(lm) {}
    ~SphUSP() {}

    double compute_coeff(const int lx, const int ly) const {
      mpfr::mpreal::set_default_prec(GMPPREC);
      const int am = abs(angular_momentum_[1]);
      const int j = (lx + ly - am) / 2;
      if ((lx + ly - am) % 2 != 0 || j < 0) {
        return 0.0;
      } else {
        const mpreal zero = "0.0";
        const mpreal pi = static_cast<mpreal>(atan(1) * 4);
        const int lz = angular_momentum_[0] - lx - ly;
        Factorial fact;
        const int lmam = angular_momentum_[0] - am;
        const mpreal lmamf = fact(lmam);
        const mpreal lpamf = fact(angular_momentum_[0] + am);
        const mpreal lf = fact(angular_momentum_[0]);
        const mpreal prefactor = sqrt(0.5 * (2*angular_momentum_[0]+1) * lmamf / lpamf / pi) / std::pow(2, angular_momentum_[0]) / lf;
        const int parity = (am - lx) % 2;
        double factor;
        if (angular_momentum_[1] > 0 && parity == 0) {
          factor = 1.0;
        } else if (angular_momentum_[1] == 0 && lx % 2 == 0) {
          factor = 1.0 / std::sqrt(2.0);
        } else if (angular_momentum_[1] < 0 && parity != 0) {
          factor = 1.0;
        } else {
          factor = 0.0;
        }

        mpreal si = zero;
        const int jp = lmam / 2;
        for (int i = j; i <= jp; ++i) {
          si += comb(angular_momentum_[0], i) * comb(i, j) * std::pow(-1.0, i)
              * ((2 * angular_momentum_[0] - 2 * i) >= 0 ? fact(2 * angular_momentum_[0] - 2 * i) : zero)
              / ((angular_momentum_[0] - am - 2 * i) >= 0 ? fact(angular_momentum_[0] - am - 2 * i) : zero);
        }
        mpreal sk = zero;
        for (int k = 0; k <= j; ++k) {
          if (lx - 2*k >= 0 && lx - 2*k <= am) {
            sk += (k <= j ? comb(j, k) : zero) * comb(am, lx - 2*k) * std::pow(-1.0, static_cast<int>((am - lx + 2*k) / 2));
          }
        }
        return (prefactor * si * sk * factor).toDouble();
      }
    }

    void expand() {
      int cnt = 0;
      for (int lz = 0; lz <= angular_momentum_[0]; ++lz) {
        for (int ly = 0; ly <= angular_momentum_[0] - lz; ++ly) {
          ++cnt;
          const int lx = angular_momentum_[0] - lz - ly;
          const double coeff = compute_coeff(lx, ly);
          std::pair<double, int> c_usp(coeff, cnt);
          usp_.push_back(c_usp);
        }
      }
    }

    int angular_momentum(const int i) const { return angular_momentum_[i]; }

    std::vector<std::pair<double, int>> usp() const { return usp_; }

    void print() {
      std::cout << "** Real spherical harmonics to unitary sphere polynomials **" << std::endl;
      std::cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                               << angular_momentum_[1] << ")" << std::endl;
      std::cout << "Y_lm = sum_i c_i * usp_i" << std::endl;
      this->expand();
      for (auto& it : usp_) {
        std::cout << "(" << std::setw(17) << std::setprecision(9) << it.first << ", " << it.second << ")" << std::endl;
      }
    }
};

class CarSph {
  protected:
    int angular_momentum_;
    std::array<int, 3> exp_;
    std::vector<double> spherical_;
    std::vector<double> cartesian_;

  public:
    CarSph(const std::array<int, 3> exp) : exp_(exp) {
      angular_momentum_ = exp_[0] + exp_[1] + exp_[2];
      for (int i = 0; i != 2 * angular_momentum_ + 1; ++i) {
        spherical_.push_back(0.0);
      }
    }
    ~CarSph() {}


    void transform_CarSph() {
      const int LEND = 7;
      const int carsphindex = LEND * angular_momentum_;
      const static bagel::CarSphList carsphlist;
      const int dimcar = (angular_momentum_ + 2) * (angular_momentum_ + 1) / 2;

      int index = 0;
      for (int i = 0; i != exp_[0]; ++i) index += angular_momentum_ + 1 - i;
      index += exp_[1];
      for (int i = 0; i != dimcar; ++i) {
        if (i - index == 0) {
          cartesian_.push_back(1.0);
        } else {
          cartesian_.push_back(0.0);
        }
      }

      double* car = cartesian_.data();
      double* sph = spherical_.data();
      carsphlist.carsphfunc_call(carsphindex, 1, car, sph);
    }

    std::vector<double> spherical() const { return spherical_; }

    std::vector<double> cartesian() const { return cartesian_; }

    int angular_momentum() const { return angular_momentum_; }

    void print() {
      std::cout << "-- Print Cartesian in carsph --" << std::endl;
      std::cout << "car[" << (angular_momentum_ + 2) * (angular_momentum_ + 1) / 2 << "]" << std::endl;
      for (auto it = cartesian_.begin(); it != cartesian_.end(); ++it) {
        std::cout << std::setw(17) << std::setprecision(9) << *it << std::endl;
      }
      std::cout << "-- Print Spherical in carsph --" << std::endl;
      std::cout << "sph[" << 2 * angular_momentum_ + 1 << "]" << std::endl;
      int cnt = 0;
      for (auto it = spherical_.begin(); it != spherical_.end(); ++it) {
        if (cnt % 2 == 0) {
          std::cout << "(" << angular_momentum_ << ", " << std::setw(3) << angular_momentum_ - cnt/2 << ") " << std::setw(17) << std::setprecision(9) << *it << std::endl;
        } else {
          std::cout << "(" << angular_momentum_ << ", " << std::setw(3) << -(angular_momentum_ - cnt/2) << ") " << std::setw(17) << std::setprecision(9) << *it << std::endl;
        }
        ++cnt;
      }
    }

};

class SH {
  protected:

  public:

  SH() {}
  ~SH() {}

  double alegendre(const int l, const int m, const double x) const {
    if (m < 0 || m > l || fabs(x) > 1.0) throw std::runtime_error("SH: m must be in [0, l] and x in [-1, 1]");
    double pmm = 1.0;
    if (m > 0) {
      double somx2 = std::sqrt((1.0 - x)*(1.0 + x));
      double fact = 1.0;
      for (int i = 1; i <= m; ++i) {
        pmm *= -fact * somx2;
        fact += 2.0;
      }
    }
    if (l == m) {
      return pmm;
    } else {
      double pmmp1 = x * (2.0 * m + 1) * pmm;
      if (l == m+1) {
        return pmmp1;
      } else {
        double plm = 0.0;
        for (int i = m + 2; i <= l; ++i) {
          plm = (x * (2 * i -1) * pmmp1 - (i + m - 1) * pmm) / (i - m);
          pmm = pmmp1;
          pmmp1 = plm;
        }
        return plm;
      }
    }
  }

  std::complex<double> ylm(const int l, const int m, const double theta, const double phi) const {
    const mpreal pi = static_cast<mpreal>(atan(1) * 4);
    Factorial fact;
    const double cth = cos(theta);
    if (abs(m) > l) throw std::runtime_error ("SH.ylm: |m| > l");

    const int am = abs(m);
    const double plm = alegendre(l, am, cth);
    mpreal f = fact(l-am)/fact(l+am);
    const double coef = std::sqrt((2*l+1) * f.toDouble() * 0.25/pi.toDouble());
    double real = coef * plm * cos(am * phi);
    double imag = coef * plm * sin(am * phi);
    if (m < 0) {
      real *= std::pow(-1, m);
      imag *= std::pow(-1, m+1);
    }
    std::complex<double> y (real, imag);
    return y;
  }

  double realSH(const int l, const int m, const double theta, const double phi) const {
    const mpreal pi = static_cast<mpreal>(atan(1) * 4);
    Factorial fact;
    const double cth = cos(theta);
    if (abs(m) > l) throw std::runtime_error ("SH.ylm: |m| > l");

    const int am = abs(m);
    const double plm = alegendre(l, am, cth);
    mpreal f = fact(l-am) / fact(l+am);
    const double coef = std::sqrt((2*l+1) * f.toDouble() * 0.25/pi.toDouble());
    if (m == 0) {
      return coef * plm;
    } else if (m > 0) {
      return std::sqrt(2.0) * std::pow(-1.0, m) * coef * plm * cos(am*phi);
    } else {
      return std::sqrt(2.0) * std::pow(-1.0, m) * coef * plm * sin(am*phi);
    }
  }

};

class RealSH : public SH {
  protected:

    std::array<int, 2> angular_momentum_;
    std::array<double, 3> centre_;

  public:

    explicit RealSH(std::array<int, 2> lm, std::array<double, 3> c) :
      angular_momentum_(lm),
      centre_(c)
      {}
    ~ RealSH() {}

    double centre(const int i) const { return centre_[i]; }

    int angular_momentum(const int i) const { return angular_momentum_[i]; }

    void print() {
      std::cout << "** Real Spherical Harmonics **" << std::endl;
      std::cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                               << angular_momentum_[1] << ")" << std::endl;
      std::cout << "Centre = " ;
      for (int i = 0; i != 3; ++i) {
        std::cout << std::setw(17) << std::setprecision(9) << centre_[i] << ";   ";
      }
      std::cout << std::endl;
      std::cout << "---" << std::endl;
    }

};

class CartesianGauss {
  protected:

    double norm_;
    double exponent_;
    double weight_;
    std::array<int, 3> angular_momentum_;
    std::array<double, 3> centre_;

  public:

    explicit CartesianGauss(const double alpha, const std::array<int, 3> l, const std::array<double, 3> c) :
      exponent_(alpha),
      angular_momentum_(l),
      centre_(c)
      { norm_ = normalise(); }
    ~CartesianGauss() {}

    double normalise() const {
      Factorial fact;
      const mpreal pi = static_cast<mpreal>(atan(1) * 4);
      const mpreal mexp = static_cast<mpreal>(exponent_);
      const int l = angular_momentum_[0] + angular_momentum_[1] + angular_momentum_[2];
      const mpreal fact1 = fact(angular_momentum_[0]) * fact(angular_momentum_[1]) * fact(angular_momentum_[2]);
      const mpreal fact2 = fact(2*angular_momentum_[0]) * fact(2*angular_momentum_[1]) * fact(2*angular_momentum_[2]);
      const mpreal mnorm = static_cast<mpreal>(pow(2 * mexp / pi, 0.75) * sqrt(pow(4 * mexp, l) * pow(2, l) * fact1 / fact2));
      return mnorm.toDouble();
    }

    int angular_momentum(const int i) const { return angular_momentum_[i]; }

    double centre(const int i) const { return centre_[i]; }

    double exponent () const { return exponent_; }

    void set_weight(const double wt) { weight_ = wt; }

    double compute(const std::array<double, 3> centre) {
      std::array<double, 3> RA;
      for (int i = 0; i != 3; ++i) RA[i] = centre[i] - centre_[i];
      double rsq = RA[0] * RA[0] + RA[1] * RA[1] + RA[2] * RA[2];
      return norm_ * std::pow(RA[0], angular_momentum_[0]) * std::pow(RA[1], angular_momentum_[1]) * std::pow(RA[2], angular_momentum_[2]) * std::exp(-exponent_ * rsq);
    }

    double compute1D(const int iX, const double centreX) {
      Factorial fact;
      const mpreal pi = static_cast<mpreal>(atan(1) * 4);
      const mpreal mexp = static_cast<mpreal>(exponent_);
      const int lx = angular_momentum_[iX];
      const mpreal mnormX = static_cast<mpreal>(pow(2 * mexp / pi, 0.25) * sqrt(pow(4 * mexp, lx) * pow(2, lx) * fact(lx) / fact(2*lx)));
      double r = centreX - centre_[iX];
      return mnormX.toDouble() * std::pow(r, lx) * std::exp(-exponent_ * r * r);
    }

    void print() {
      std::cout << "** Cartesian Gaussian **" << std::endl;
      std::cout << "Angular momentum (lx, ly, lz) = (" << angular_momentum_[0] << ", "
                                                       << angular_momentum_[1] << ", "
                                                       << angular_momentum_[2] << ")" << std::endl;
      std::cout << "Centre = " ;
      for (int i = 0; i != 3; ++i) {
        std::cout << std::setw(17) << std::setprecision(9) << centre_[i] << ";   ";
      }
      std::cout << std::endl;
      std::cout << "Exponent alpha = " << std::setw(17) << std::setprecision(9) << exponent_  << std::endl;
      std::cout << "Normalisation const N = " << std::setw(17) << std::setprecision(9) << norm_ << std::endl;
      std::cout << "---" << std::endl;
    }

};

class ProjectionInt {
  protected:

    std::shared_ptr<const CartesianGauss> gauss_;
    std::shared_ptr<const RealSH> sh_;

  public:

    ProjectionInt(std::shared_ptr<const CartesianGauss> gauss, std::shared_ptr<const RealSH> sh) : gauss_(gauss), sh_(sh) {}
    ~ProjectionInt() {}

    double integrate3SHs(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3) const {
      const double pi = static_cast<double>(atan(1.0) * 4.0);
      static Wigner3j wigner3j;
      const double w1 = wigner3j.lookup_wigner3j(l1, m1, l2, m2, l3, m3);
      const double w2 = wigner3j.lookup_wigner3j(l1, 0, l2, 0, l3, 0);
      const double coeff = std::sqrt(0.25 * (2*l1 + 1) * (2*l2 + 1) * (2*l3 + 1) / pi);

      return coeff * w1 * w2;
    }

    double integrate3USP(const int i, const int j, const int k) const {
      const mpreal pi = static_cast<mpreal>(atan(1.0) * 4.0);
      if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
        Factorial fact;
        mpreal num = "2.0" * fact(i) * fact(j) * fact(k) * fact((i+j+k+2)/2);
        mpreal denom = fact(i/2) * fact(j/2) * fact(k/2) * fact(i+j+k+2);
        return ("4.0" * pi * num / denom).toDouble();
      } else {
        return 0.0;
      }
    }

    double integrate2SH1USP(const std::pair<int, int> lm1, const std::pair<int, int> lm2, const std::array<int, 3> ijk) const {
      std::vector<std::pair<double, int>> usp;

      std::array<int, 2> lm = {lm1.first, lm1.second};
      std::shared_ptr<SphUSP> sphusp = std::make_shared<SphUSP>(lm);
      int cnt = 0;
      int nonzero = 0;
      for (int lz = 0; lz <= lm1.first; ++lz) {
        for (int ly = 0; ly <= lm1.first - lz; ++ly) {
          const int lx = lm1.first - lz - ly;
          const double coeff = sphusp->compute_coeff(lx, ly);
          if (coeff != 0.0) {
            nonzero ++;
            std::pair<double, int> c_usp(coeff, cnt);
            usp.push_back(c_usp);
          }
          ++cnt;
        }
      }
      const int n1 = nonzero;

      lm = {lm2.first, lm2.second};
      sphusp = std::make_shared<SphUSP>(lm);
      cnt = 0;
      nonzero = 0;
      for (int lz = 0; lz <= lm2.first; ++lz) {
        for (int ly = 0; ly <= lm2.first - lz; ++ly) {
          const int lx = lm2.first - lz - ly;
          const double coeff = sphusp->compute_coeff(lx, ly);
          if (coeff != 0.0) {
            nonzero ++;
            std::pair<double, int> c_usp(coeff, cnt);
            usp.push_back(c_usp);
          }
          ++cnt;
        }
      }
      const int n2 = nonzero;

      assert (n1 + n2 == usp.size());
      double ans = 0.0;
      for (int i = 0; i != n1; ++i) {
        for (int j = n1; j != n1 + n2; ++j) {
          const double coeff = usp[i].first * usp[j].first;
          std::array<int, 3> ki;
          int id = usp[i].second;
          int kz = 0;
          for (int lp1 = lm1.first + 1; lp1 != 0; --lp1) {
            if (id - lp1 < 0) {
              ki[2] = kz;
              ki[1] = id;
              ki[0] = lm1.first - ki[2] - ki[1];
              break;
            } else {
              kz++;
              id -= lp1;
            }
          }
          std::array<int, 3> kj;
          id = usp[j].second;
          kz = 0;
          for (int lp1 = lm2.first + 1; lp1 != 0; --lp1) {
            if (id - lp1 < 0) {
              kj[2] = kz;
              kj[1] = id;
              kj[0] = lm2.first - kj[2] - kj[1];
              break;
            } else {
              kz++;
              id -= lp1;
            }
          }
          const int x = ki[0] + kj[0] + ijk[0];
          const int y = ki[1] + kj[1] + ijk[1];
          const int z = ki[2] + kj[2] + ijk[2];
          ans += coeff * integrate3USP(x, y, z);
        }
      }

      return ans;

    }

    double compute(const double r) const {
      const double pi = static_cast<double>(atan(1.0) * 4.0);
      std::array<double, 3> AB;
      for (int i = 0; i != 3; ++i) AB[i] = gauss_->centre(i) - sh_->centre(i);
        const double dAB = std::sqrt(AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]);
        const double thAB = acos(AB[2]/dAB);
        const double phAB = atan2(AB[1], AB[0]);
        const int nx = gauss_->angular_momentum(0);
        const int ny = gauss_->angular_momentum(1);
        const int nz = gauss_->angular_momentum(2);
        const int nu = nx + ny + nz;
        const int lnu = sh_->angular_momentum(0) + nu;
        const mpreal exponential = static_cast<mpreal>(exp(-gauss_->exponent() * (dAB - r) * (dAB - r)));
        double ans = 0.0;
        for (int kx = 0; kx <= nx; ++kx) {
          const double ckx = comb(nx, kx) * std::pow(AB[0], nx - kx);
          for (int ky = 0; ky <= ny; ++ky) {
            const double cky = comb(ny, ky) * std::pow(AB[1], ny - ky);
            for (int kz = 0; kz <= nz; ++kz) {
              const double ckz = comb(nz, kz) * std::pow(AB[2], nz - kz);
              const int lk = kx + ky + kz;
              const double rxyz = std::pow(r, lk);

              double sld = 0.0;
              for (int ld = 0; ld <= lnu; ++ld) {
                double smu = 0.0;
                for (int m = 0; m <= 2 * ld; ++m) {
                  const int mu = m - ld;
                  const double Z_AB = (dAB == 0 ? (1.0/std::sqrt(4.0*pi)) : sh_->realSH(ld, mu, thAB, phAB));

                  const std::array<int, 3> exp = {kx, ky, kz};
                  const std::pair<int, int> lm1(ld, mu);
                  const std::pair<int, int> lm2(sh_->angular_momentum(0), sh_->angular_momentum(1));
                  smu += Z_AB * integrate2SH1USP(lm1, lm2, exp);
//                std::cout << "ld = " << ld << " mu = " << mu << " thAB = " << thAB << " phAB = " <<  phAB << " Z_AB = " << Z_AB << " 2SH1USP = " << integrate2SH1USP(lm1, lm2, exp) << std::endl;
                }
                Modified_Spherical_Bessel_Iexp msbessel(ld);
                const mpreal sbessel = msbessel.compute(2.0 * gauss_->exponent() * dAB * r);
                sld += smu * sbessel.toDouble();
              }
              ans += sld * ckx * cky * ckz * rxyz * std::pow(-1.0, lk - nu);
            }
          }
        }
//      return ans * 4.0 * pi * gauss_->normalise() * exponential.toDouble();
        return ans * 4.0 * pi * exponential.toDouble();
    }

};

class ECP_Type2 {

  protected:
    std::shared_ptr<ProjectionInt> projAB_;
    std::shared_ptr<ProjectionInt> projCB_;
    int nkl_; // 0, -1, -2
    double zetakl_;

  public:

    ECP_Type2(std::tuple<int, double, std::shared_ptr<ProjectionInt>, std::shared_ptr<ProjectionInt>> ang) :
      nkl_(std::get<0>(ang)),
      zetakl_(std::get<1>(ang)),
      projAB_(std::get<2>(ang)),
      projCB_(std::get<3>(ang))
    {}

    ~ECP_Type2() {}

    mpreal compute(const mpreal r) {
      const double projABr = projAB_->compute(r.toDouble());
      const double projCBr = projCB_->compute(r.toDouble());
      return static_cast<mpreal>(projABr * pow(r, nkl_ + 2) * exp(-zetakl_ * r * r) * projCBr * r * r);
    }

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
      const mpreal pi = static_cast<mpreal>(atan(1) * 4);
      mpreal mbessi = "0.0";
      Factorial fact;
      const mpreal mfact = fact(l);
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

}

#endif
