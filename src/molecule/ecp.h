//
// BAGEL - Parallel electron correlation program.
// Filename: ecp.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_MOLECULE_ECP_H
#define __SRC_MOLECULE_ECP_H

#include <vector>
#include <src/molecule/shell_ECP.h>

namespace bagel {

class ECP {

  protected:
    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp_;

  public:
    ECP(std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp) : shells_ecp_(shells_ecp) {}
    ~ECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp() const { return shells_ecp_; }
    std::shared_ptr<const Shell_ECP> shell_ecp(const int i) const { return shells_ecp_[i]; }

    std::shared_ptr<const Shell_ECP> shell_maxl_ecp() const {
      int maxl = 0;
      std::shared_ptr<const Shell_ECP> maxl_shell;
      for (auto& ish : shells_ecp_) {
        if (ish->angular_number() >= maxl) {
          maxl = ish->angular_number();
          maxl_shell = ish;
        }
      }

      return maxl_shell;

    }

    void print() const {
      for (auto& i : shells_ecp_) std::cout << i->show() << std::endl;
    }

};

}

#endif

