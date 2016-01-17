//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ecp.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __SRC_MOLECULE_ECP_H
#define __SRC_MOLECULE_ECP_H

#include <vector>
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

  private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & ecp_ncore_ & ecp_maxl_ & shells_ecp_ & ishell_maxl_ & nshell_ & nr_;
    }

  public:
    ECP();
    ECP(const int ncore, const int maxl, std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp);
    ~ECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_ecp() const { return shells_ecp_; }
    std::shared_ptr<const Shell_ECP> shell_ecp(const int i) const { return shells_ecp_[i]; }

    void get_shell_maxl_ecp();

    std::shared_ptr<const Shell_ECP> shell_maxl_ecp() const;
    int ecp_maxl() const { return ecp_maxl_; }

    int ecp_ncore() const { return ecp_ncore_; }

    int nshell() const { return nshell_; }

    const std::array<int, 3> nr() { return nr_; }
    int nr(const int i) const { return nr_[i]; }

    double position(const int i) const { return shells_ecp_[0]->position(i); };
    const std::array<double,3>& position() const { return shells_ecp_[0]->position(); };

    void print() const;

};

}

#endif

