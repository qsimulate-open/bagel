//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fmminfo.h
// Copyright (C) 2017 Toru Shiozaki
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


#ifndef __SRC_WFN_FMMINFO_H
#define __SRC_WFN_FMMINFO_H

#include <src/molecule/shellpair.h>
#include <src/molecule/atom.h>

namespace bagel {

class FMMInfo {

  protected:
    std::vector<std::shared_ptr<const ShellPair>> shellpairs_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & shellpairs_;
    }

  public:
    FMMInfo() { }
    FMMInfo(const std::vector<std::shared_ptr<const Atom>>& atoms,
            const std::vector<std::vector<int>>& offsets, const std::string extent_type = "yang");
    ~FMMInfo() { }

    std::vector<std::shared_ptr<const ShellPair>> shellpairs() const { return shellpairs_; }
    std::shared_ptr<const ShellPair> shellpair(const int i) const { return shellpairs_[i]; }
    int nshellpair() const { return shellpairs_.size(); }
};

}

#endif
