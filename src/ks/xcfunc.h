//
// BAGEL - Parallel electron correlation program.
// Filename: xcfunc.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_KS_XCFUNC_H
#define __SRC_KS_XCFUNC_H

#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <config.h>
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
      map_.insert(std::make_pair("b3lyp", XC_HYB_GGA_XC_B3LYP));
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
      std::cout << "      * DFT functional " << name_ << " using libxc" << std::endl << std::endl;
    }
    ~XCFunc() { xc_func_end(&func_); }

    std::unique_ptr<double[]> compute_exc(int np, const std::unique_ptr<double[]>& rho, const std::unique_ptr<double[]>& sigma) const {
      std::unique_ptr<double[]> out(new double[np]);
      if (func_.info->family == XC_FAMILY_HYB_GGA || func_.info->family == XC_FAMILY_GGA) { 
        xc_gga_exc(&func_, np, rho.get(), sigma.get(), out.get());
      } else {
        throw std::runtime_error("So far only GGA and Hybrid GGA is supported.");
      }
      return std::move(out);
    }

};

#else
class XCFunc { XCFunc(const std::string) { assert(false); } }; // dummy
#endif

}

#endif
