//
// BAGEL - Parallel electron correlation program.
// Filename: vertex_sp.h
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


#ifndef __SRC_PERIODIC_VERTEX_SP_H
#define __SRC_PERIODIC_VERTEX_SP_H

#include <src/util/constants.h>
#include <src/molecule/shellpair.h>

namespace bagel {

class VertexSP {
  protected:
    std::bitset<64> key_;
    std::shared_ptr<const ShellPair> sp_;
    std::vector<std::shared_ptr<const ZMatrix>> multipole_;

  public:
    VertexSP(std::bitset<64> key, std::shared_ptr<const ShellPair> sp, const int lmax = 10)
     : key_(key), sp_(sp) { multipole_ = sp_->multipoles(lmax); }
    ~VertexSP() { }

    std::bitset<64> key() const { return key_; }
    std::bitset<3> node_key(const int i) const {
      std::bitset<3> out;
      out[0] = key_[i * 3    ];
      out[1] = key_[i * 3 + 1];
      out[2] = key_[i * 3 + 2];
      return out;
    }

    std::array<double, 3> centre() const { return sp_->centre(); }
    double centre(const int i) const { return sp_->centre(i); }
    std::shared_ptr<const ShellPair> sp() const { return sp_; }
    std::array<std::shared_ptr<const Shell>, 2> shells() const { return sp_->shells(); }
    std::shared_ptr<const Shell> shell0() const { return sp_->shell(0); }
    std::shared_ptr<const Shell> shell1() const { return sp_->shell(1); }
    bool is_neighbour(std::shared_ptr<const VertexSP> v, const int ws) const { return sp_->is_neighbour(v->sp(), ws); }
    int shell_ind(const int i) const { assert(i==0 || i==1); return (i == 0) ? sp_->shell_ind(0) : sp_->shell_ind(1); }
    std::array<int, 2> offset() const { return sp_->offset(); }
    int offset(const int i) const { assert(i==0 || i==1); return sp_->offset(i); }
    double extent() const { return sp_->extent(); }
    int nbasis0() const { return sp_->nbasis0(); }
    int nbasis1() const { return sp_->nbasis1(); }
    std::vector<std::shared_ptr<const ZMatrix>> multipole() const { return multipole_; }
    int nmult() const { return multipole_.size(); }
    double schwarz() const { return sp_->schwarz(); }
};

}
#endif
