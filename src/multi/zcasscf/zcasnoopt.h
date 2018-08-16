//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasnoopt.h
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#ifndef __SRC_MULTI_ZCASSCF_ZCASNOOPT_H
#define __SRC_MULTI_ZCASSCF_ZCASNOOPT_H

#include <src/multi/zcasscf/zcasscf.h>

namespace bagel {

class ZCASNoopt_base : public ZCASSCF {
  protected:
    virtual void init_mat1e() override final { /*do nothing*/ }
    virtual void impose_symmetry(std::shared_ptr<ZMatrix>) const override final { }
    virtual void impose_symmetry(std::shared_ptr<ZRotFile>) const override final { }
    virtual bool kramers() const = 0;

  protected:
    ZCASNoopt_base(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : ZCASSCF(idat, geom, ref) { }

  public:
    void compute() override final;
};


class ZCASNoopt : public ZCASNoopt_base {
  protected:
    virtual void init_coeff() override;
    virtual bool kramers() const override { return true; }

  public:
    ZCASNoopt(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref);
    std::shared_ptr<const Reference> conv_to_ref() const override { return conv_to_ref_(true); }
};


class ZCASNoopt_London : public ZCASNoopt_base {
  protected:
    virtual void init_coeff() override;
    virtual bool kramers() const override { return false; }

  public:
    ZCASNoopt_London(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref);
    std::shared_ptr<const Reference> conv_to_ref() const override { return conv_to_ref_(false); }
};

}

#endif
