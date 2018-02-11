//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: jop_london.h
// Copyright (C) 2017 Toru Shiozaki
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

#ifndef __BAGEL_CI_ZFCI_JOP_LONDON_H
#define __BAGEL_CI_ZFCI_JOP_LONDON_H

#include <src/ci/zfci/zmofile.h>

namespace bagel {

class Jop_London : public ZMOFile {
  protected:
    std::shared_ptr<ZMatrix> compute_hcore() const override;
    std::shared_ptr<ZMatrix> compute_fock(std::shared_ptr<const ZMatrix> hcore, const int nclosed, const bool, const bool) const override;
    std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;
    std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;

  public:
    Jop_London(const std::shared_ptr<const Geometry> geom, const int nstart, const int nfence, std::shared_ptr<const ZCoeff_Block> coeff);
};

}

#endif
