//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/apply_operator.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_RAS_APPLY_OPERATOR_H
#define __SRC_RAS_APPLY_OPERATOR_H

#include <src/ci/ciutil/citraits.h>
#include <src/ci/ras/civector.h>
#include <src/asd/dmrg/product_civec.h>
#include <src/asd/gamma_sq.h>

namespace bagel {

/// helper class used to apply second quantization operators used in plays like ASD_DMRG and similar in spirit to FormSigmaRAS
class ApplyOperator {
  private:
    /// takes a potential target string, returns whether it is valid and flips the proper bit in the process
    bool test_and_set(std::bitset<nbit__>& bit, const bool action, const int orbital) {
      bool out = action ? bit[orbital] : !bit[orbital];
      action ? bit.reset(orbital) : bit.set(orbital);
      return out;
    }

    std::shared_ptr<const RASBlock<double>> get_block(const RASCivecView& source, const std::array<int, 4>& dhp,
                                                      const std::shared_ptr<const RASBlock<double>>& target) const;
    std::shared_ptr<const CIBlockInfo<RASString>> get_blockinfo(const std::shared_ptr<const RASDeterminants>& sourcedet, const std::array<int, 4>& dhp,
                                                          const std::shared_ptr<const CIBlockInfo<RASString>>& target) const;

  public:
    ApplyOperator() {}

    /// \f[ |\mbox{target}\rangle \leftarrow a \hat E_{\mbox{operations}} |\mbox{source}\rangle + |\mbox{target}\rangle \f]
    void operator()(const double a, const RASCivecView source, RASCivecView target, const std::vector<GammaSQ>& operations, const std::vector<int>& orbitals) const;

    /// \f[ |\mbox{target}\rangle \leftarrow a \hat E_{\mbox{operations}} |\mbox{source}\rangle + |\mbox{target}\rangle \f]
    void operator()(const double a, const RASBlockVectors& source, RASBlockVectors& target, const std::vector<GammaSQ>& operations, const std::vector<int>& orbitals) const;
};

}

#endif
