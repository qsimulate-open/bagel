//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_MP2_MP2_H
#define __SRC_MP2_MP2_H

#include <src/scf/hf/rhf.h>
#include <src/wfn/method.h>

namespace bagel {

class MP2 : public Method {
  protected:
    std::shared_ptr<RHF> scf_;
    int ncore_;

    std::string abasis_;

    double energy_;

  public:
    MP2(const std::shared_ptr<const PTree>, const std::shared_ptr<const Geometry>, const std::shared_ptr<const Reference> = nullptr);

    virtual void compute() override;
    virtual std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    double energy() const { return energy_; }
    int ncore() const { return ncore_; }
    std::string abasis() const { return abasis_; }
    std::shared_ptr<const RHF> scf() const { return scf_; }
};

}

#endif
