//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison.h
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

// Desc :: The implementation closely follows Harrison and Zarrabian
//

#ifndef __BAGEL_ZFCI_ZHARRISONZARRABIAN_H
#define __BAGEL_ZFCI_ZHARRISONZARRABIAN_H

#include <src/zfci/zfci.h>

namespace bagel {

class ZHarrison : public ZFCI {

  protected:
    void const_denom() override { }

    std::shared_ptr<RelMOFile> jop_; // this hides ZFCI::jop_

    // run-time functions.
#if 0
    void sigma_aa(std::shared_ptr<const ZDvec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const MOFile> jop) const { }
    void sigma_bb(std::shared_ptr<const ZDvec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const MOFile> jop) const { }
    void sigma_2ab_1(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const { }
    void sigma_2ab_2(std::shared_ptr<ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const MOFile> jop) const { }
    void sigma_2ab_3(std::shared_ptr<ZCivec> sigma, std::shared_ptr<ZDvec> e) const { }
#endif

  public:
    // this constructor is ugly... to be fixed some day...
    ZHarrison(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
        const int ncore = -1, const int nocc = -1, const int nstate = -1) : ZFCI(a, g, b, ncore, nocc, nstate) {
      relupdate();
    }

    std::shared_ptr<RelZDvec> form_sigma(std::shared_ptr<const RelZDvec> c, std::shared_ptr<const ZMOFile_base> jop, const std::vector<int>& conv) const override {
      return std::shared_ptr<RelZDvec>();
    }

    void update(std::shared_ptr<const Coeff>) override { assert(false); }

    void relupdate() {
      Timer timer;
      jop_ = std::make_shared<RelJop>(ref_, ncore_, ncore_+norb_*2, "HZ");

      // right now full basis is used.
      std::cout << "    * Integral transformation done. Elapsed time: " << std::setprecision(2) << timer.tick() << std::endl << std::endl;
      const_denom();
    }
};

}

#endif

