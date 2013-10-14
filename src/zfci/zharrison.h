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

#ifndef __BAGEL_ZFCI_ZHARRISON_H
#define __BAGEL_ZFCI_ZHARRISON_H

#include <src/wfn/method.h>
#include <src/zfci/relmofile.h>
#include <src/zfci/reldvec.h>
#include <src/rel/relreference.h>

namespace bagel {

class ZHarrison : public Method {

  protected:
    // max #iteration
    int max_iter_;
    // threshold for variants
    double thresh_;
    double print_thresh_;

    // numbers of electrons
    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;

    // number of states
    int nstate_;

    // total energy
    std::vector<double> energy_;

    // CI vector at convergence
    std::shared_ptr<RelZDvec> cc_;

    // MO integrals
    std::shared_ptr<RelMOFile> jop_;

    // Determinant space
    std::shared_ptr<const RelSpace> space_;
    std::shared_ptr<const RelSpace> int_space_;

    // denominator
    std::shared_ptr<RelDvec> denom_;

    // some init functions
    void common_init(); // may end up unnecessary
    // obtain determinants for guess generation
    void generate_guess(const int nspin, const int nstate, std::shared_ptr<RelZDvec>);
    // generate spin-adapted guess configurations
    virtual std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet);

    // print functions
    void print_header() const;

    void const_denom();

    // run-time functions.
    // aaaa and bbbb
    void sigma_aa(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const RelMOFile> jop, const bool trans = false) const;
    void sigma_1e_ab(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const RelMOFile> jop, const bool trans = false) const;

    void sigma_2e_annih_ab(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;
    void sigma_2e_annih_aa(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;

    void sigma_2e_create_bb(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;
    void sigma_2e_create_ab(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;

    void sigma_2e_h0101_h1001(std::shared_ptr<const ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const RelMOFile> jop) const;

    template<int i, int j, int k, int l>
    void sigma_2e_h(std::shared_ptr<const ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const RelMOFile> jop, const bool trans) const {
      static_assert(!((i|1)^1) && !((j|1)^1) && !((k|1)^1) && !((l|1)^1), "illegal call of sigma_2e_h");
      const int ij = d->ij();
      const int lenab = d->lena()*d->lenb();
      std::stringstream ss; ss << i << j << k << l;
      std::bitset<4> bit4(ss.str());
      if (trans) bit4 = ~bit4;
      zgemm3m_("n", "t", lenab, ij, ij, 1.0, d->data(), lenab, jop->mo2e(bit4)->data(), ij, 0.0, e->data(), lenab);
    }

    void sigma_one(std::shared_ptr<const ZCivec> cc, std::shared_ptr<RelZDvec> sigmavec, std::shared_ptr<const RelMOFile> jop,
                   const int istate, const bool diag, const bool trans) const;

  public:
    // this constructor is ugly... to be fixed some day...
    ZHarrison(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
        const int ncore = -1, const int nocc = -1, const int nstate = -1);

    std::shared_ptr<RelZDvec> form_sigma(std::shared_ptr<const RelZDvec> c, std::shared_ptr<const RelMOFile> jop, const std::vector<int>& conv) const;

    void update() {
      Timer timer;
      jop_ = std::make_shared<RelJop>(ref_, ncore_, ncore_+norb_*2);

      // right now full basis is used.
      std::cout << "    * Integral transformation done. Elapsed time: " << std::setprecision(2) << timer.tick() << std::endl << std::endl;
      const_denom();
    }

    void compute() override;

    // returns members
    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    double core_energy() const { return jop_->core_energy(); }

    int nij() const { return norb_*norb_; }

    // TODO
    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }
};

}

#endif

