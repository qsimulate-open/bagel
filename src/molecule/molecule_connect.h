//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molecule_connect.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __BAGEL_MOLECULE_MOLECULE_CONNECT_H
#define __BAGEL_MOLECULE_MOLECULE_CONNECT_H

#include <set>
#include <src/molecule/molecule.h>

namespace bagel {
namespace molecule_details {

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


bool is_one_molecule(std::list<std::shared_ptr<Node>> nodes)  {
  std::list<std::shared_ptr<Node>> disconnected = nodes;
  std::list<std::shared_ptr<Node>> connected;
  std::list<std::shared_ptr<Node>> newly_added;
  connected.push_back(*nodes.begin());
  newly_added.push_back(*nodes.begin());
  disconnected.erase(disconnected.begin());
  bool done = false;

  auto add_connected = [&connected, &disconnected, &newly_added, &done] () {
    std::list<std::shared_ptr<Node>> added;
    done = true;
    for (auto i = newly_added.begin(); i != newly_added.end(); ++i) {
      for (auto j = disconnected.begin(); j != disconnected.end(); ) {
        if ((*i)->connected_with(*j)) {
          done = false;
          added.push_back(*j);
          connected.push_back(*j);
          j = disconnected.erase(j);
        } else {
          ++j;
        }
      }
    }
    newly_added = added;
  };

  while (!done) {
    add_connected();
  }

  return (disconnected.size() == 0);
}

}
}

#endif
