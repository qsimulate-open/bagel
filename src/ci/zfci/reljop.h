//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reljop.h
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

#ifndef __BAGEL_CI_ZFCI_RELJOP_H
#define __BAGEL_CI_ZFCI_RELJOP_H

#include <src/ci/zfci/zmofile.h>
#include <src/df/reldffull.h>

namespace bagel {

class RelJop : public ZMOFile {
  protected:
    bool gaunt_;
    bool breit_;

    virtual std::shared_ptr<ZMatrix> compute_hcore() const override;
    std::shared_ptr<ZMatrix> compute_fock(std::shared_ptr<const ZMatrix> hcore, const int nclosed, const bool store_c, const bool store_g) const override;
    std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;
    std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;

  public:
    RelJop(const std::shared_ptr<const Geometry> geom, const int nstart, const int nfence, std::shared_ptr<const ZCoeff_Block> coeff,
           const bool gaunt, const bool breit, const bool store_c = false, const bool store_g = false);

    static std::tuple<std::list<std::shared_ptr<RelDFHalf>>, std::list<std::shared_ptr<RelDFHalf>>>
      compute_half(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit);
    static std::shared_ptr<ListRelDFFull> compute_full(std::shared_ptr<const ZMatrix> coeff, std::list<std::shared_ptr<RelDFHalf>> half, const bool appj);
    static std::shared_ptr<ListRelDFFull> compute_full(std::shared_ptr<const ZMatrix> coeff, std::list<std::shared_ptr<const RelDFHalf>> half, const bool appj);
};

}

#endif
