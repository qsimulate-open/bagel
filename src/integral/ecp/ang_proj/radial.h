//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#ifndef __SRC_INTEGRAL_ECP_ANG_PROJ_RADIAL_H
#define __SRC_INTEGRAL_ECP_ANG_PROJ_RADIAL_H

#include "proj.h"

using namespace mpfr;

template<typename T, typename... Value>
class Radial_Int {

  protected:
    int ngrid_;
    int max_iter_;
    double thresh_int_;
    vector<mpreal> x_;
    vector<mpreal> w_;
    T function_;

  public:
    Radial_Int(const int max_iter, const double thresh_int, Value... tail) :
      max_iter_(max_iter),
      thresh_int_(thresh_int),
      function_(tail...)
    {
      integrate();
    }

    ~Radial_Int() {}

    void integrate() {
      double ans = 0.0;
      double previous = 0.0;
      int ngrid = 31;
      for (int iter = 0; iter != max_iter_; ++iter) {
        GaussChebyshev2nd(ngrid);
        vector<mpreal> r;
        transform_Becke(r);
        ans = 0.0;
        int cnt = 0;
        for (auto& it : r) {
          ans += (function_.compute(it) * w_[cnt]).toDouble();
          ++cnt;
        }
        const double error = ans - previous;
        std::cout << "Iteration no. " << iter << " ngrid = " << ngrid << " ans = " << ans << " error = " << error << std::endl;
        if (error < thresh_int_ && iter != 0) {
          std::cout << "Integration converged..." << std::endl;
          std::cout << "Radial integral = " << ans << std::endl;
          break;
        } else if (iter == max_iter_-1) {
          std::cout << "Max iteration exceeded..." << std::endl;
        }
        previous = ans;
        x_.clear();
        w_.clear();
        ngrid *= 2;
      }
    }

    void transform_Log3(vector<mpreal>& r) {

    }

    void transform_Becke(vector<mpreal>& r) {
      const double alpha = 1.8;
      const mpreal one = "1.0";
      for (auto& it : x_) {
        r.push_back(static_cast<mpreal>(alpha * (one + it) / (one - it)));
      }
    }

    void GaussChebyshev2nd(const int ngrid) {
      const mpreal pi = static_cast<mpreal>(atan(1) * 4);
      for (int i = 1; i != ngrid; ++i) {
        x_.push_back(static_cast<mpreal>(cos(i*pi/(ngrid+1))));
        w_.push_back(static_cast<mpreal>(pi * sin(i*pi/(ngrid+1)) * sin(i*pi/(ngrid+1)) / (ngrid+1)));
      }
    }

};

#endif
