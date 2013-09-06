//
// BAGEL - Parallel electron correlation program.
// Filename: geometry_connect.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __BAGEL_WFN_GEOMETRY_CONNECT_H
#define __BAGEL_WFN_GEOMETRY_CONNECT_H

#include <set>
#include <src/wfn/geometry.h>

namespace bagel {
namespace geometry_details {

class Node {
  protected:
    const std::shared_ptr<const Atom> myself_;
    const int num_;
    // in order to avoid cyclic references
    std::list<std::weak_ptr<Node>> connected_;


  public:
    Node(const std::shared_ptr<const Atom> o, const int n) : myself_(o), num_(n) {};

    void add_connected(const std::shared_ptr<Node> i) {
      std::weak_ptr<Node> in = i;
      for (auto& iter : connected_)
        if (iter.lock() == i) throw std::logic_error("Node::add_connected");
      connected_.push_back(in);
    }

    const std::shared_ptr<const Atom> atom() const { return myself_; }
    int num() const { return num_; }

    std::set<std::shared_ptr<Node>> common_center(const std::shared_ptr<Node> o) const {
      std::set<std::shared_ptr<Node>> out;
      for (auto& c : connected_) {
        if (c.lock()->connected_with(o)) out.insert(c.lock());
      }
      return out;
    }

    bool connected_with(const std::shared_ptr<Node> o) {
      bool out = false;
      for (auto& i : connected_) {
        if (i.lock() == o) {
          out = true;
          break;
        }
      }
      return out;
    }

};

static double adf_rho(std::shared_ptr<const Node> i, std::shared_ptr<const Node> j) {
  return std::exp(- i->atom()->distance(j->atom()) / (i->atom()->cov_radius()+j->atom()->cov_radius()) + 1.0);
}

}
}

#endif
