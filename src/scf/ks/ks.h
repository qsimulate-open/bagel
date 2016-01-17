//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ks.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#ifndef __BAGEL_SRC_KS_KS_H
#define __BAGEL_SRC_KS_KS_H

#include <src/scf/scf_base.h>
#include <src/scf/ks/xcfunc.h>
#include <src/scf/ks/dftgrid.h>

// I only implement a DF version

namespace bagel {

class KS : public SCF_base {
  protected:
    std::string name_;
    std::shared_ptr<XCFunc> func_;
    std::shared_ptr<DFTGrid_base> grid_;

  public:
    KS(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re = nullptr)
      : SCF_base(idata, geom, re) {

      std::cout << indent << "*** Kohn-Sham DFT ***" << std::endl << std::endl;

      // default is now B3LYP
      name_ = idata->get<std::string>("xc_func", "b3lyp");
      func_ = std::make_shared<XCFunc>(name_);

      Timer preptime;
      grid_ = std::make_shared<DefaultGrid>(geom);
      preptime.tick_print("DFT grid generation");

      std::cout << std::endl;

      if (re) std::cout << " **** we have not implemented DFT with a reference ****";
    }

    virtual void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    std::shared_ptr<const XCFunc> func() const { return func_; }
    std::shared_ptr<const DFTGrid_base> grid() const { return grid_; }
};

}

#endif
