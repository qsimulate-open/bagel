//
// BAGEL - Parallel electron correlation program.
// Filename: smith_info.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __SRC_SMITH_SMITH_INFO_H
#define __SRC_SMITH_SMITH_INFO_H

#include <src/wfn/reference.h>

namespace bagel {

class SMITH_Info : public Reference {
  protected:
    int ncore_;
    double thresh_;
    int maxiter_;

  public:
    SMITH_Info(std::shared_ptr<const Reference> o, const std::shared_ptr<const PTree> idata) : Reference(*o) {
      const bool frozen = idata->get<bool>("frozen", true);
      ncore_ = idata->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
      if (ncore_)
        std::cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << std::endl;

      thresh_ = idata->get<double>("thresh", 1.0e-8);
      maxiter_ = idata->get<int>("maxiter", 50);
    }

    int ncore() const { return ncore_; }
    double thresh() const { return thresh_; }
    int maxiter() const { return maxiter_; }
};

}

#endif
