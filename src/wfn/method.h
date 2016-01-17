//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: method.h
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

#ifndef __SRC_WFN_METHOD_H
#define __SRC_WFN_METHOD_H

#include <src/wfn/geometry.h>
#include <src/wfn/reference.h>

// this file should be header only (in order not to introduce additional dependency)

namespace bagel {

class Method {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & idata_ & geom_ & ref_;
    }

  public:
    Method() { }
    Method(std::shared_ptr<const PTree> p, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r)
     : idata_(p), geom_(g), ref_(r) { }
    virtual ~Method() { }

    virtual void compute() = 0;
    virtual std::shared_ptr<const Reference> conv_to_ref() const = 0;

    std::shared_ptr<const PTree> idata() const { return idata_; }
    std::shared_ptr<const Reference> ref() const { return ref_; }
    std::shared_ptr<const Geometry> geom() const { return geom_; }

};

}

#endif
