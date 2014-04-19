//
// BAGEL - Parallel electron correlation program.
// Filename: radial.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_INTEGRAL_ECP_RADIAL_H
#define __SRC_INTEGRAL_ECP_RADIAL_H

#include <cmath>
#include <vector>
#include <src/util/constants.h>

namespace bagel {

template<typename T, typename... Value>
class RadialInt {

  protected:
    int ngrid_;
    bool print_intermediate_;
    int max_iter_;
    double thresh_int_;
    std::vector<double> x_;
    std::vector<double> w_;
    double integral_;

  public:
    RadialInt(Value... tail, const bool print = false, const int max_iter = 100, const double thresh_int = 1e-10) :
      print_intermediate_(print),
      max_iter_(max_iter),
      thresh_int_(thresh_int)
    {
      T function(tail...);
      integrate(function);
    }

    ~RadialInt() {}

    void integrate(T function) {
      double ans = 0.0;
      double previous = 0.0;
      int ngrid = 31;
      for (int iter = 0; iter != max_iter_; ++iter) {
        GaussChebyshev2nd(ngrid);
        std::vector<double> r;
        transform_Becke(r);
        ans = 0.0;
        int cnt = 0;
        for (auto& it : r) {
          ans += function.compute(it) * w_[cnt];
          ++cnt;
        }
        const double error = ans - previous;
        if (print_intermediate_)
           std::cout << "Iteration no. " << iter << " ngrid = " << ngrid << " ans = " << ans << " error = " << error << std::endl;
        if (fabs(error) < thresh_int_ && iter != 0) {
          if (print_intermediate_) {
            std::cout << "Integration converged..." << std::endl;
            std::cout << "Radial integral = " << ans << std::endl;
          }
          integral_ = ans;
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

    double integral() { return integral_; }

    void transform_Becke(std::vector<double>& r) {
      const double alpha = 1.0;
      int cnt = 0;
      for (auto& it : x_) {
        r.push_back(alpha * (1.0 + it) / (1.0 - it));
        w_[cnt] *= 2.0 / (1.0 - it) / (1.0 - it);
        ++cnt;
      }
    }

    void GaussChebyshev2nd(const int ngrid) {
      for (int i = 1; i != ngrid; ++i) {
        x_.push_back(cos(i * pi__ / (ngrid + 1)));
        w_.push_back(pi__ * sin(i * pi__ / (ngrid + 1)) / (ngrid + 1));
      }
    }

};

}

#endif
