//
// BAGEL - Parallel electron correlation program.
// Filename: branch.h
// Copyright (C) 2016 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_BRANCH_H
#define __SRC_PERIODIC_BRANCH_H

#include <src/molecule/shellpair.h>

namespace bagel {

class Branch {
  friend class Box;
  protected:
    int iws_;
    std::array<double, 3> centre_;
    std::vector<std::shared_ptr<const ShellPair>> sp_;

    double extent_;
    std::vector<std::shared_ptr<const ShellPair>> neigh_;
    std::vector<std::shared_ptr<const ShellPair>> non_neigh_;

    bool is_neigh(std::shared_ptr<const Branch> b) const;
    void get_neigh(const std::vector<std::shared_ptr<Branch>>& branch);

  public:
    Branch(const int ws, const std::array<double, 3>& c, const std::vector<std::shared_ptr<const ShellPair>>& sp);
    ~Branch() { }

    double extent() const { return extent_; }
    int ws() const { return iws_; }
    const std::array<double, 3>& centre() const { return centre_; }
    const std::vector<std::shared_ptr<const ShellPair>>& neigh() const { return neigh_; }
    const std::vector<std::shared_ptr<const ShellPair>>& non_neigh() const { return non_neigh_; }
    int nsp() const { return sp_.size(); }
    const std::vector<std::shared_ptr<const ShellPair>>& sp() const { return sp_; }
    std::shared_ptr<const ShellPair> sp(const int i) const { return sp_[i]; }
};

}
#endif
