//
// BAGEL - Parallel electron correlation program.
// Filename: relmofile.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell  <caldwell@u.northwestern.edu>
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


#ifndef __BAGEL_ZFCI_RELMOFILE_H
#define __BAGEL_ZFCI_RELMOFILE_H

#include <unordered_map>
#include <src/math/zmatrix.h>
#include <src/rel/reldffull.h>
#include <src/rel/dfock.h>
#include <src/rel/relreference.h>

namespace bagel {

class RelMOFile {
  protected:
    int nocc_;
    int nbasis_;
    double core_energy_;

    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const ZMatrix> core_fock_;
    std::shared_ptr<const ZMatrix> coeff_;
    std::array<std::shared_ptr<const ZMatrix>,2> kramers_coeff_;

    bool gaunt_;
    bool breit_;

    // creates integral files and returns the core energy.
    void init(const int nstart, const int nend);

    // hamiltoniam data
    std::unordered_map<std::bitset<2>, std::shared_ptr<const ZMatrix>> mo1e_;
    std::unordered_map<std::bitset<4>, std::shared_ptr<const ZMatrix>> mo2e_;

    // generates Kramers symmetry-adapted orbitals
    void compress_and_set(std::unordered_map<std::bitset<2>, std::shared_ptr<const ZMatrix>> buf1e,
                          std::unordered_map<std::bitset<4>, std::shared_ptr<const ZMatrix>> buf2e);

    virtual std::unordered_map<std::bitset<2>, std::shared_ptr<const ZMatrix>> compute_mo1e(const std::array<std::shared_ptr<const ZMatrix>,2> coeff) = 0;
    virtual std::unordered_map<std::bitset<4>, std::shared_ptr<const ZMatrix>> compute_mo2e(const std::array<std::shared_ptr<const ZMatrix>,2> coeff) = 0;

    // half transformed integrals for CASSCF
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_coulomb_;
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_gaunt_;

  public:
    RelMOFile(const std::shared_ptr<const Geometry>, std::shared_ptr<const ZMatrix>, const bool gaunt, const bool breit);

    // static function
    static std::array<std::shared_ptr<const ZMatrix>,2> kramers(std::shared_ptr<const ZMatrix> coeff, std::shared_ptr<const ZMatrix> overlap, std::shared_ptr<const ZMatrix> eig);


    std::shared_ptr<const ZMatrix> core_fock() const { return core_fock_; }

    std::shared_ptr<const ZMatrix> mo1e(const std::bitset<2>& b) const { return mo1e_.at(b); }
    std::shared_ptr<const ZMatrix> mo2e(const std::bitset<4>& b) const { return mo2e_.at(b); }
//  const std::complex<double>& mo1e(const std::bitset<2>& b, const size_t i) const { return mo1e_.at(b)->data(i); }
    const std::complex<double>& mo1e(const std::bitset<2>& b, const size_t i, const size_t j) const { return mo1e_.at(b)->element(i,j); }
    const std::complex<double>& mo2e(const std::bitset<4>& b, const size_t i, const size_t j, const size_t k, const size_t l) const { return mo2e_.at(b)->element(i+nocc_*j, k+nocc_*l); }
    std::shared_ptr<const ZMatrix> mo1e(std::string&& b) const { return mo1e(std::bitset<2>(std::move(b))); }
    std::shared_ptr<const ZMatrix> mo2e(std::string&& b) const { return mo2e(std::bitset<4>(std::move(b))); }
    const std::complex<double>& mo1e(std::string&& b, const size_t i, const size_t j) const { return mo1e(std::bitset<2>(std::move(b)), i, j); }
    const std::complex<double>& mo2e(std::string&& b, const size_t i, const size_t j, const size_t k, const size_t l) const { return mo2e(std::bitset<4>(std::move(b)), i, j, k, l); }

    double core_energy() const { return core_energy_; }

    std::array<std::shared_ptr<const ZMatrix>,2> kramers_coeff() const { return kramers_coeff_; }
    std::shared_ptr<const ZMatrix> coeff() const {
      auto coeff_tot = std::make_shared<ZMatrix>(kramers_coeff_[0]->ndim(), nocc_*2);
      coeff_tot->copy_block(0,     0, kramers_coeff_[0]->ndim(), nocc_, kramers_coeff_[0]);
      coeff_tot->copy_block(0, nocc_, kramers_coeff_[1]->ndim(), nocc_, kramers_coeff_[1]);
      return coeff_tot;
    }
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_coulomb() const { return half_complex_coulomb_; }
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_gaunt() const { return half_complex_gaunt_; }


    static std::unordered_map<std::bitset<2>, std::shared_ptr<const RelDFFull>>
        compute_full(std::array<std::array<std::shared_ptr<const Matrix>,4>,2> rocoeff, std::array<std::array<std::shared_ptr<const Matrix>,4>,2> iocoeff,
                     std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half, const bool appj, const bool appjj = false) {
      std::unordered_map<std::bitset<2>, std::shared_ptr<const RelDFFull>> out;
      for (size_t t = 0; t != 4; ++t) {
        std::list<std::shared_ptr<RelDFFull>> dffull;
        for (auto& i : half[t/2])
          dffull.push_back(std::make_shared<RelDFFull>(i, rocoeff[t%2], iocoeff[t%2]));
        DFock::factorize(dffull);
        assert(dffull.size() == 1);
        dffull.front()->scale(dffull.front()->fac()); // take care of the factor
        out[std::bitset<2>(t)] = dffull.front();
        assert(!appj || !appjj);
        if (appj)
          out.at(std::bitset<2>(t)) = out.at(std::bitset<2>(t))->apply_J();
        else if (appjj)
          out.at(std::bitset<2>(t)) = out.at(std::bitset<2>(t))->apply_JJ();
      }
      return out;
    }
};


class RelJop : public RelMOFile {
  protected:
    std::unordered_map<std::bitset<2>, std::shared_ptr<const ZMatrix>> compute_mo1e(const std::array<std::shared_ptr<const ZMatrix>,2> coeff) override;
    std::unordered_map<std::bitset<4>, std::shared_ptr<const ZMatrix>> compute_mo2e(const std::array<std::shared_ptr<const ZMatrix>,2> coeff) override;

  public:
    RelJop(const std::shared_ptr<const Geometry> geo, const int c, const int d, std::shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit)
      : RelMOFile(geo, coeff, gaunt, breit) { init(c, d); }
};


#if 0
class RelHtilde : public ZHtilde_Base, public RelMOFile {
  protected:
    std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) override { return std::make_tuple(h1_tmp_, 0.0); };
    std::shared_ptr<const ZMatrix> compute_mo2e(const int, const int) override { return h2_tmp_; };
  public:
    RelHtilde(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const ZMatrix> h1, std::shared_ptr<const ZMatrix> h2)
      : ZHtilde_Base(h1, std::move(h2)), RelMOFile(b) {
      core_energy_ = create_Jiiii(c, d);
    }
};
#endif

}

#endif
