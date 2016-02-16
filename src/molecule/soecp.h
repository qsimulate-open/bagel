//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: soecp.h
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

  private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & shells_so_ & so_maxl_;
    }

  public:
    SOECP() : shells_so_(1, std::make_shared<const Shell_ECP>()), so_maxl_(0) {}
    SOECP(std::vector<std::shared_ptr<const Shell_ECP>> shells_so) : shells_so_(shells_so) {
#ifndef NDEBUG
      for (auto& i : shells_so)
        assert(i->angular_number() > 0);
#endif
      so_maxl_ = shells_so.back()->angular_number();
    }
    ~SOECP() {}

    std::vector<std::shared_ptr<const Shell_ECP>> shells_so() const { return shells_so_; }
    std::shared_ptr<const Shell_ECP> shell_so(const int i) const { return shells_so_[i]; }

    int nshell() const { return shells_so_.size(); }
    int so_maxl() const { return so_maxl_; }

    double position(const int i) const { return shells_so_.front()->position(i); };
    const std::array<double,3>& position() const { return shells_so_.front()->position(); };

    void print() const {
      std::cout << "+++ SOECP Parameters +++" << std::endl;
      for (auto& i : shells_so_) std::cout << i->show() << std::endl;
    }
};

}

#endif

