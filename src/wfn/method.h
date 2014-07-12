//
// BAGEL - Parallel electron correlation program.
// Filename: method.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_WFN_METHOD_H
#define __SRC_WFN_METHOD_H

#include <src/wfn/geometry.h>
#include <src/wfn/geometry_london.h>
#include <src/wfn/reference.h>

// this file should be header only (in order not to introduce additional dependency)

namespace bagel {

class Method_ {
  protected:
    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Reference> ref_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & idata_ & ref_;
    }

  public:
    Method_() { }
    Method_(std::shared_ptr<const PTree> p, std::shared_ptr<const Reference> r)
     : idata_(p), ref_(r) { }
    virtual ~Method_() { }

    virtual void compute() = 0;
    virtual std::shared_ptr<const Reference> conv_to_ref() const = 0;

    std::shared_ptr<const PTree> idata() const { return idata_; }
    std::shared_ptr<const Reference> ref() const { return ref_; }

};


template<typename GeomType>
class SubMethod : public Method_ {
  protected:
    std::shared_ptr<const GeomType> geom_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Method_);
      ar & geom_;
    }

  public:
    SubMethod() { }
    SubMethod(std::shared_ptr<const PTree> p, std::shared_ptr<const GeomType> g, std::shared_ptr<const Reference> r)
     : Method_(p, r), geom_(g) { }
    virtual ~SubMethod() { }

    std::shared_ptr<const GeomType> geom() const { return geom_; }
};

using Method = SubMethod<Geometry>;
using Method_London = SubMethod<Geometry_London>;

}

#endif
