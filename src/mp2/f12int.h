//
// BAGEL - Parallel electron correlation program.
// Filename: f12int.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_MP2_F12INT_H
#define __SRC_MP2_F12INT_H

#include <memory>
#include <map>
#include <src/mp2/f12mat.h>
#include <src/wfn/geometry.h>
#include <src/wfn/reference.h>

namespace bagel {

class F12Int {
  protected:

    const std::multimap<std::string, std::string> idata_;
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const Reference> ref_;
    double gamma_;

    std::shared_ptr<F12Mat> bmat_;
    std::shared_ptr<F12Mat> xmat_;
    std::shared_ptr<F12Mat> vmat_;

  public:
    F12Int(const std::multimap<std::string, std::string>, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref,
           const double, const int);
    ~F12Int() {};

    std::shared_ptr<F12Mat> robust_fitting(std::shared_ptr<const DFFullDist> doo, std::shared_ptr<const DFFullDist> yoo);

};

}


#endif

