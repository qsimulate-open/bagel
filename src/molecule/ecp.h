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
#include <algorithm>
#include <src/molecule/shell_ECP.h>

namespace bagel {

class ECP {

  protected:
    int ecp_ncore_;
    int ecp_maxl_;
    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp_;
    int ishell_maxl_;
    int nshell_;
    std::array<int, 3> nr_;

  public:
    ECP(const int ncore, const int maxl, std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp)
     : ecp_ncore_(ncore), ecp_maxl_(maxl), shells_ecp_(shells_ecp) {
      get_shell_maxl_ecp();
      nshell_ = shells_ecp_.size();
    }
    ~ECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp() const { return shells_ecp_; }
    std::shared_ptr<const Shell_ECP> shell_ecp(const int i) const { return shells_ecp_[i]; }

    void get_shell_maxl_ecp() {
      for (auto ish = shells_ecp_.begin(); ish != shells_ecp_.end(); ++ish)
        if ((*ish)->angular_number() == ecp_maxl_) {
          ishell_maxl_ = std::distance(shells_ecp_.begin(), ish);
          break;
        }

      for (int i = 0; i != 3; ++i)
        nr_[i] = std::count(shells_ecp_[ishell_maxl_]->ecp_r_power().begin(), shells_ecp_[ishell_maxl_]->ecp_r_power().end(), std::abs(i-2));

    }

    std::shared_ptr<const Shell_ECP> shell_maxl_ecp() const { return shells_ecp_[ishell_maxl_]; }
    const int ecp_maxl() const { return ecp_maxl_; }

    const int ecp_ncore() const { return ecp_ncore_; }

    const int nshell() const { return nshell_; }

    const std::array<int, 3> nr() { return nr_; }
    const int nr(const int i) const { return nr_[i]; }

    void print() const {
      std::cout << "+++ ECP Parameters +++" << std::endl;
      std::cout << "Number of core electrons = " << ecp_ncore_ << std::endl;
      for (auto& i : shells_ecp_) std::cout << i->show() << std::endl;
    }

};

}

#endif

