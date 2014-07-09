//
// BAGEL - Parallel electron correlation program.
// Filename: ecp.h
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


#ifndef __SRC_MOLECULE_ECP_H
#define __SRC_MOLECULE_ECP_H

#include <vector>
#include <memory>
#include <algorithm>
#include <src/molecule/shellecp.h>

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
    ECP();
    ECP(const int ncore, const int maxl, std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp);
    ~ECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp() const { return shells_ecp_; }
    std::shared_ptr<const Shell_ECP> shell_ecp(const int i) const { return shells_ecp_[i]; }

    void get_shell_maxl_ecp();

    std::shared_ptr<const Shell_ECP> shell_maxl_ecp() const;
    const int ecp_maxl() const { return ecp_maxl_; }

    const int ecp_ncore() const { return ecp_ncore_; }

    const int nshell() const { return nshell_; }

    const std::array<int, 3> nr() { return nr_; }
    const int nr(const int i) const { return nr_[i]; }

    double position(const int i) const { return shells_ecp_[0]->position(i); };
    const std::array<double,3>& position() const { return shells_ecp_[0]->position(); };

    void print() const;

};

}

#endif

