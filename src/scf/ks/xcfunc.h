//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: xcfunc.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SRC_KS_XCFUNC_H
#define __SRC_KS_XCFUNC_H

#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <bagel_config.h>
#ifdef HAVE_XC_H
#include <xc.h> // header provided by libxc
#endif

namespace bagel {

#ifdef HAVE_XC_H
class FuncList {
  protected:
    std::map<std::string, int> map_;
  public:
    FuncList() {
      map_.emplace("slater",    XC_LDA_X);
      map_.emplace("xalpha",    XC_LDA_C_XALPHA);
      map_.emplace("pw92c",     XC_LDA_C_PW);
      map_.emplace("pbex",      XC_GGA_X_PBE);
      map_.emplace("pbec",      XC_GGA_C_PBE);
      map_.emplace("b3lyp",     XC_HYB_GGA_XC_B3LYP);
      map_.emplace("pbe0",      XC_HYB_GGA_XC_PBEH);
    }
    int num(const std::string& name) const {
      auto iter = map_.find(name);
      if (iter == map_.end())
        throw std::runtime_error("Functional  " + name + " not found");
      return iter->second;
    }
};


class XCFunc {
  protected:
    xc_func_type func_;
    const std::string name_;
    const FuncList flist;

  public:
    XCFunc(const std::string name) : name_(name) {
      if (xc_func_init(&func_, flist.num(name_), XC_UNPOLARIZED))
        throw std::runtime_error("unknown functional..");
    }
    ~XCFunc() { xc_func_end(&func_); }

    void compute_exc_vxc(int np, const double* rho, const double* sigma, double* exc, double* vxc, double* vxc2) const {
      if (lda()) {
        xc_lda_exc_vxc(&func_, np, rho, exc, vxc);
      } else if (gga()) {
        xc_gga_exc_vxc(&func_, np, rho, sigma, exc, vxc, vxc2);
      } else {
        throw std::runtime_error("Meta GGA not supported yet");
      }
    }

    void compute_vxc(int np, const double* rho, const double* sigma, double* vxc, double* vxc2) const {
      if (lda()) {
        xc_lda_vxc(&func_, np, rho, vxc);
      } else if (gga()) {
        xc_gga_vxc(&func_, np, rho, sigma, vxc, vxc2);
      } else {
        throw std::runtime_error("Meta GGA not supported yet");
      }
    }

    bool lda() const { return func_.info->family == XC_FAMILY_LDA; }
    bool gga() const { return func_.info->family == XC_FAMILY_HYB_GGA || func_.info->family == XC_FAMILY_GGA; }

    double scale_ex() const { return (func_.info->family == XC_FAMILY_HYB_GGA) ? xc_hyb_exx_coef(&func_) : 0.0; }
};

#else
class XCFunc {
public:
  XCFunc(const std::string) { assert(false); }
  void compute_exc_vxc(int np, const double* rho, const double* sigma, double* exc, double* vxc, double* vxc2) const {}
  void compute_vxc(int np, const double* rho, const double* sigma, double* vxc, double* vxc2) const {}
  bool lda() const { return true; }
  double scale_ex() const { return 0.0; }
}; // dummy
#endif

}

#endif
