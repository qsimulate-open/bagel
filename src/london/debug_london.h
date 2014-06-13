//
// BAGEL - Parallel electron correlation program.
// Filename: debug_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#if 0
#ifndef __BAGEL_SRC_LONDON_DEBUG_LONDON_H
#define __BAGEL_SRC_LONDON_DEBUG_LONDON_H

#include <src/wfn/method.h>

namespace bagel {

class Debug_London : public Method {

  protected:

  public:
    Debug_London() { }
    Debug_London(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry_London> cgeom, const std::shared_ptr<const Reference> re = nullptr);
    virtual ~Debug_London() { }

    void compute() override;
    std::shared_ptr<const Reference> conv_to_ref() const override;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Debug_London)

#endif
#endif
