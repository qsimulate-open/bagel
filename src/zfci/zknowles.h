//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
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

// Desc :: The implementation closely follows Knowles and Handy 1984 CPL.
//         It is amazing how easy it is to implement FCI !!
//

#ifndef __BAGEL_ZFCI_ZKNOWLESHANDY_H
#define __BAGEL_ZFCI_ZKNOWLESHANDY_H

#include <src/zfci/zfci.h>

namespace bagel {

class ZKnowlesHandy : public ZFCI {
  protected:
    // denominator
    void const_denom() override;

    // virtual application of Hamiltonian
    std::shared_ptr<RelZDvec> form_sigma(std::shared_ptr<const RelZDvec> c, std::shared_ptr<const ZMOFile_Base> jop, const std::vector<int>& conv) const override;

    // run-time functions
    void sigma_1(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZMOFile_Base> jop) const;
    void sigma_3(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZMOFile_Base> jop) const;
    void sigma_2b (std::shared_ptr<ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const ZMOFile_Base> jop) const;
    void sigma_2c1(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;
    void sigma_2c2(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;

    // multiply copmlex phase factors to each orbitals (TODO only for debugging)
    void mult_phase_factor();

  public:
    // this constructor is ugly... to be fixed some day...
    ZKnowlesHandy(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
                  const int ncore = -1, const int nocc = -1, const int nstate = -1) : ZFCI(a, g, b, ncore, nocc, nstate) {
      auto relref = std::dynamic_pointer_cast<const RelReference>(ref_);
      if (!relref) {
        update(ref_->coeff());
      } else {
        relupdate();
      }
    }
    void relupdate();
    void update(std::shared_ptr<const Coeff>) override;
};

}

#endif

