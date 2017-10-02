//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relmofile.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell  <caldwell@u.northwestern.edu>
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


#ifndef __BAGEL_ZFCI_RELMOFILE_H
#define __BAGEL_ZFCI_RELMOFILE_H

#include <unordered_map>
#include <src/util/kramers.h>
#include <src/util/math/zmatrix.h>
#include <src/df/reldffull.h>
#include <src/wfn/relreference.h>

namespace bagel {

class RelMOFile {
  protected:
    int nocc_;
    int nbasis_;
    double core_energy_;

    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const ZMatrix> core_fock_;
    std::shared_ptr<const RelCoeff_Block> coeff_;
    std::shared_ptr<Kramers<1,ZMatrix>> kramers_coeff_;

    bool gaunt_;
    bool breit_;

    // creates integral files and returns the core energy.
    void init(const int nstart, const int nfence, const bool store_c, const bool store_g);

    // hamiltoniam data
    std::shared_ptr<Kramers<2,ZMatrix>> mo1e_;
    std::shared_ptr<Kramers<4,ZMatrix>> mo2e_;

    // generates Kramers symmetry-adapted orbitals
    void compress_and_set(std::shared_ptr<Kramers<2,ZMatrix>> buf1e,
                          std::shared_ptr<Kramers<4,ZMatrix>> buf2e);

    virtual std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) = 0;
    virtual std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) = 0;

  public:
    RelMOFile(const std::shared_ptr<const Geometry>, std::shared_ptr<const RelCoeff_Block>,
              const bool gaunt, const bool breit);

    std::shared_ptr<const ZMatrix> core_fock() const { return core_fock_; }

    template<typename T>
    std::shared_ptr<const ZMatrix> mo1e(const T& b) const { KTag<2> bb(b); return mo1e_->at(bb); }
    template<typename T>
    std::shared_ptr<const ZMatrix> mo2e(const T& b) const { KTag<4> bb(b); return mo2e_->at(bb); }
    template<typename T>
    const std::complex<double>& mo1e(const T& b, const size_t i, const size_t j) const { return mo1e(b)->element(i,j); }
    template<typename T>
    const std::complex<double>& mo2e(const T& b, const size_t i, const size_t j, const size_t k, const size_t l) const { return mo2e(b)->element(i+nocc_*j, k+nocc_*l); }

    std::shared_ptr<const Kramers<2,ZMatrix>> mo1e() const { return mo1e_; }
    std::shared_ptr<const Kramers<4,ZMatrix>> mo2e() const { return mo2e_; }

    double core_energy() const { return core_energy_; }

    std::shared_ptr<const RelCoeff_Block> coeff() const { return coeff_; }

    static std::tuple<std::list<std::shared_ptr<RelDFHalf>>, std::list<std::shared_ptr<RelDFHalf>>>
      compute_half(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit);
    static std::shared_ptr<ListRelDFFull> compute_full(std::shared_ptr<const ZMatrix> coeff, std::list<std::shared_ptr<RelDFHalf>> half, const bool appj);
    static std::shared_ptr<ListRelDFFull> compute_full(std::shared_ptr<const ZMatrix> coeff, std::list<std::shared_ptr<const RelDFHalf>> half, const bool appj);
};


class RelJop : public RelMOFile {
  protected:
    std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;
    std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<1,ZMatrix>> coeff) override;

  public:
    RelJop(const std::shared_ptr<const Geometry> geom, const int nstart, const int nfence, std::shared_ptr<const RelCoeff_Block> coeff,
      const bool gaunt, const bool breit, const bool store_c = false, const bool store_g = false)
      : RelMOFile(geom, coeff, gaunt, breit) { init(nstart, nfence, store_c, store_g); }
};

}

#endif
