//
// BAGEL - Parallel electron correlation program.
// Filename: shellpair.h
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


#ifndef __SRC_MOLECULE_SHELLPAIR_H
#define __SRC_MOLECULE_SHELLPAIR_H

#include <src/molecule/shell.h>

namespace bagel {

class ShellPair {

  protected:
    std::array<std::shared_ptr<const Shell>, 2> shells_;
    std::array<int, 2> offset_;
    std::pair<int, int> shell_ind_;
    std::string extent_type_;
    int nbasis0_, nbasis1_;

    double thresh_;
    double schwarz_;
    std::array<double, 3> centre_;
    double extent_;
    void init();

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & shells_ & offset_ & shell_ind_ & extent_type_ & nbasis0_ & nbasis1_
         & thresh_ & schwarz_ & centre_ & extent_;
    }

  public:
    ShellPair() { }
    ShellPair(const std::array<std::shared_ptr<const Shell>, 2>& shells, const std::array<int, 2>& offset,
              const std::pair<int, int>& shell_ind, const std::string extent_type = "yang", const double thresh = 1e-10);
    bool is_neighbour(std::shared_ptr<const ShellPair> sp, const double ws) const;

    const std::array<std::shared_ptr<const Shell>, 2>& shells() const { return shells_; }
    std::shared_ptr<const Shell> shell(const int i) const { assert(i==0 || i==1); return shells_[i]; }
    const std::array<int, 2>& offset() const { return offset_; }
    int offset(const int i) const { assert(i==0 || i==1); return offset_[i]; }
    const std::pair<int, int>& shell_ind() const { return shell_ind_; }
    int shell_ind(const int i) const { assert(i==0 || i==1); return (i == 0) ? shell_ind_.first : shell_ind_.second; }
    double schwarz() const { return schwarz_; }
    const std::array<double, 3>& centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }
    double extent() const { return extent_; }
    int nbasis0() const { return nbasis0_; }
    int nbasis1() const { return nbasis1_; }

    std::vector<std::shared_ptr<const ZMatrix>> multipoles(int lmax = 10, const std::array<double, 3>& Qcentre = {{0,0,0}}) const;
};

}

#endif
