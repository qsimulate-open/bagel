//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldffullt.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_DF_RELDFFULLT_H
#define __SRC_DF_RELDFFULLT_H

#include <src/df/dfdistt.h>
#include <src/df/reldffull.h>

namespace bagel {

class RelDFFullT {
  protected:
    std::shared_ptr<const SpinorInfo> basis_;
    std::array<std::shared_ptr<DFDistT>,2> dffull_;

  public:
    RelDFFullT(std::shared_ptr<const RelDFFull> full, std::shared_ptr<const StaticDist> dist = nullptr);

    const double* data(const int i) const { return dffull_.at(i)->data(); }

    // to make this class compatible with DFDistT
    int naux()   const { assert(dffull_[0]->naux() == dffull_[1]->naux());     return dffull_[0]->naux(); }
    int bstart() const { assert(dffull_[0]->bstart() == dffull_[1]->bstart()); return dffull_[0]->bstart(); }
    void discard_df();

    std::vector<std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>> get_slice(const int start, const int end) const;
    std::shared_ptr<ZMatrix> replicate() const;

    int locate(const size_t i, const size_t n) const { assert(dffull_[0]->locate(i,n) == dffull_[1]->locate(i,n)); return dffull_[0]->locate(i,n); }

    std::shared_ptr<const SpinorInfo> basis() const { return basis_; }

    static int nblocks() { return 1; }
};


class ListRelDFFullT {
  protected:
    std::list<std::shared_ptr<RelDFFullT>> data_;

  public:
    ListRelDFFullT(std::shared_ptr<const ListRelDFFull> full, std::shared_ptr<const StaticDist> dist = nullptr) {
      for (auto& i : *full)
        data_.push_back(std::make_shared<RelDFFullT>(i, dist));
    }

    void discard_df() {
      for (auto& i : data_)
        i->discard_df();
    }

    std::list<std::shared_ptr<RelDFFullT>>::const_iterator begin() const { return data_.cbegin(); }
    std::list<std::shared_ptr<RelDFFullT>>::const_iterator end() const { return data_.cend(); }
    std::list<std::shared_ptr<RelDFFullT>>::iterator begin() { return data_.begin(); }
    std::list<std::shared_ptr<RelDFFullT>>::iterator end() { return data_.end(); }

    std::list<std::shared_ptr<RelDFFullT>> data() { return data_; }
};

}

#endif
