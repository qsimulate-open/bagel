//
// BAGEL - Parallel electron correlation program.
// Filename: soecp.h
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


#ifndef __SRC_MOLECULE_SOECP_H
#define __SRC_MOLECULE_SOECP_H

#include <vector>
#include <algorithm>
#include <cassert>
#include <src/molecule/shellecp.h>

namespace bagel {

class SOECP {

  protected:
    std::vector<std::shared_ptr<const Shell_ECP>> shells_so_;
    int so_maxl_;

  public:
    SOECP() : shells_so_(1, std::make_shared<const Shell_ECP>()), so_maxl_(0) {}
    SOECP(std::vector<std::shared_ptr<const Shell_ECP>> shells_so) : shells_so_(shells_so) {
      for (auto& i : shells_so) assert(i->angular_number() > 0);
      so_maxl_ = shells_so.back()->angular_number();
    }
    ~SOECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_so() const { return shells_so_; }
    std::shared_ptr<const Shell_ECP> shell_so(const int i) const { return shells_so_[i]; }

    const int nshell() const { return shells_so_.size(); }
    const int so_maxl() const { return so_maxl_; }

    double position(const int i) const { return shells_so_.front()->position(i); };
    const std::array<double,3>& position() const { return shells_so_.front()->position(); };

    void print() const {
      std::cout << "+++ SOECP Parameters +++" << std::endl;
      for (auto& i : shells_so_) std::cout << i->show() << std::endl;
    }
};

}

#endif

