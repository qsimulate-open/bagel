//
// BAGEL - Parallel electron correlation program.
// Filename: soscf.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __BAGEL_SRC_SCF_SOSCF_H
#define __BAGEL_SRC_SCF_SOSCF_H

#include <src/scf/scf_base.h>
#include <src/scf/sohcore.h>

namespace bagel {

class SOSCF : public SCF_base {
  protected:
    std::shared_ptr<const Matrix> aodensity_;
    std::shared_ptr<const SOHcore> sohcore_;
    std::shared_ptr<const SOHcore_base> sohcore_base_;
    std::shared_ptr<const Coeff> socoeff_;
    std::shared_ptr<const Matrix> sooverlap_;
    std::shared_ptr<const Matrix> sotildex_;
    std::unique_ptr<double[]> soeig_;
    

  public:
    SOSCF(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry> geom,
          const std::shared_ptr<const Reference> re = std::shared_ptr<const Reference>());

    void initial_guess();

    void compute() override;

    std::shared_ptr<const Matrix> sooverlap();
    std::shared_ptr<const Matrix> sotildex();

    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }
};

}

#endif
