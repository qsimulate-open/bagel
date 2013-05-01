//
// BAGEL - Parallel electron correlation program.
// Filename: ks.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#ifndef __BAGEL_SRC_KS_KS_H
#define __BAGEL_SRC_KS_KS_H

#include <src/scf/scf_base.h>
#include <src/ks/xcfunc.h>
#include <src/ks/dftgrid.h>

// I only implement a DF version

namespace bagel {

class KS : public SCF_base {
  protected:
    std::string name_;
    std::shared_ptr<XCFunc> func_;
    std::shared_ptr<DFTGrid_base> grid_;

  public:
    KS(const boost::property_tree::ptree& idata_, const std::shared_ptr<const Geometry> geom,
        const std::shared_ptr<const Reference> re = std::shared_ptr<const Reference>())
      : SCF_base(idata_, geom, re) {

      std::cout << indent << "*** Kohn-Sham DFT ***" << std::endl << std::endl;

      // default is now B3LYP
      name_ = idata_.get<std::string>("xc_func", "b3lyp"); 
      func_ = std::make_shared<XCFunc>(name_);

      Timer preptime; 
      grid_ = std::make_shared<DefaultGrid>(geom);
      preptime.tick_print("DFT grid generation");

      std::cout << std::endl;

      if (re) throw std::runtime_error("we have not implemented DFT with a reference");
    }

    virtual void compute() override;

    std::shared_ptr<Reference> conv_to_ref() const;

    std::shared_ptr<const XCFunc> func() const { return func_; }
    std::shared_ptr<const DFTGrid_base> grid() const { return grid_; }
};

}

#endif
