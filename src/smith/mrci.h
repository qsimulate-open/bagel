//
// BAGEL - Parallel electron correlation program.
// Filename: mrci.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_SMITH_MRCI_H
#define __SRC_SMITH_MRCI_H

#include <src/wfn/reference.h>
#include <src/wfn/method.h>
#include <src/util/input/input.h>

namespace bagel {

class MRCI : public Method {
  protected:
    std::shared_ptr<const Matrix> coeff_;

    int ncore_;
    double energy_;
    double thresh_;

    std::vector<double> ref_energy_;

  public:
    MRCI(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    int ncore() const { return ncore_; }
    double energy() const { return energy_; }
    double thresh() const { return thresh_; }

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }
};

}

#endif
