//
// Newint - Parallel electron correlation program.
// Filename: mp2grad.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_MP2_MP2GRAD_H
#define __SRC_MP2_MP2GRAD_H

#include <src/mp2/mp2.h>
#include <src/wfn/reference.h>
#include <string>
#include <map>

namespace bagel {

class MP2Grad : public MP2 {
  protected:

  public:
    MP2Grad(const std::multimap<std::string, std::string>, const std::shared_ptr<const Geometry>);
    ~MP2Grad() {};

    void compute();

    std::shared_ptr<Reference> conv_to_ref() const { return ref_; };
};

}

#endif
