//
// BAGEL - Parallel electron correlation program.
// Filename: pfock.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_PFOCK_H
#define __SRC_PERIODIC_PFOCK_H

#include <src/periodic/lattice.h>
#include <src/periodic/pdata.h>

namespace bagel {

class PFock : public PData {
  protected:
    std::shared_ptr<const Lattice> lattice_;
    std::shared_ptr<const PData> previous_;
    std::shared_ptr<const PData> pcoeff_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<PData>(*this)
         & lattice_ & previous_ & pcoeff_;
    }

  public:
    PFock() { }
    PFock(std::shared_ptr<const Lattice> lattice, std::shared_ptr<const PData> previous, std::shared_ptr<const PData> pcoeff);
    ~PFock() { }

    void form_pfock();  /* add 2-electron part */
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::PFock)

#endif

