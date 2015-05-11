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
#include <src/util/kramers.h>
#include <src/util/math/zmatrix.h>
#include <src/df/reldffull.h>
#include <src/wfn/relreference.h>

namespace bagel {

class RelMOFile {
  protected:
    int nocc_;
    int nbasis_;
    int charge_;
    double core_energy_;

    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const ZMatrix> core_fock_;
    std::shared_ptr<const RelCoeff_Block> coeff_;
    std::shared_ptr<Kramers<2,ZMatrix>> kramers_coeff_;

    bool gaunt_;
    bool breit_;
    bool tsymm_;

    // creates integral files and returns the core energy.
    void init(const int nstart, const int nfence, const bool restricted = false);

    // hamiltoniam data
    std::shared_ptr<Kramers<2,ZMatrix>> mo1e_;
    std::shared_ptr<Kramers<4,ZMatrix>> mo2e_;

    // generates Kramers symmetry-adapted orbitals
    void compress_and_set(std::shared_ptr<Kramers<2,ZMatrix>> buf1e,
                          std::shared_ptr<Kramers<4,ZMatrix>> buf2e);

    virtual std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<2,ZMatrix>> coeff) = 0;
    virtual std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<2,ZMatrix>> coeff) = 0;

    // half transformed integrals for CASSCF
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_coulomb_;
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_gaunt_;

  public:
    RelMOFile(const std::shared_ptr<const Geometry>, std::shared_ptr<const RelCoeff_Block>, const int charge, const bool gaunt, const bool breit, const bool tsymm);

    // static function
    static std::shared_ptr<Kramers<2,ZMatrix>> kramers(std::shared_ptr<const ZMatrix> coeff, std::shared_ptr<const ZMatrix> overlap, std::shared_ptr<const ZMatrix> eig);

    // static function used to make the order of eigenvalues & eigenvectors of ZMatrix::diagonalize() match that given by QuatMatrix::diagonalize()
    static void rearrange_eig(VectorB& eig, std::shared_ptr<ZMatrix> coeff, const bool includes_neg = true);

    std::shared_ptr<const ZMatrix> core_fock() const { return core_fock_; }

    template<typename T>
    std::shared_ptr<const ZMatrix> mo1e(const T& b) const { KTag<2> bb(b); return mo1e_->at(bb); }
    template<typename T>
    std::shared_ptr<const ZMatrix> mo2e(const T& b) const { KTag<4> bb(b); return mo2e_->at(bb); }
    template<typename T>
    const std::complex<double>& mo1e(const T& b, const size_t i, const size_t j) const { return mo1e(b)->element(i,j); }
    template<typename T>
    const std::complex<double>& mo2e(const T& b, const size_t i, const size_t j, const size_t k, const size_t l) const { return mo2e(b)->element(i+nocc_*j, k+nocc_*l); }
    std::shared_ptr<const Kramers<4,ZMatrix>> mo2e() const { return mo2e_; }

    double core_energy() const { return core_energy_; }

    std::shared_ptr<const ZMatrix> kramers_coeff(const int i) const { return kramers_coeff_->at(i); }
    std::shared_ptr<const Kramers<2,ZMatrix>> kramers_coeff() const { return kramers_coeff_; }
    std::shared_ptr<const ZMatrix> coeff() const {
      auto coeff_tot = std::make_shared<ZMatrix>(kramers_coeff_->at(0)->ndim(), nocc_*2);
      coeff_tot->copy_block(0,     0, kramers_coeff_->at(0)->ndim(), nocc_, kramers_coeff_->at(0));
      coeff_tot->copy_block(0, nocc_, kramers_coeff_->at(1)->ndim(), nocc_, kramers_coeff_->at(1));
      return coeff_tot;
    }
    std::shared_ptr<const RelCoeff_Block> coeff_input() const { return coeff_; }
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_coulomb() const { return half_complex_coulomb_; }
    std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half_complex_gaunt() const { return half_complex_gaunt_; }


    static std::shared_ptr<const Kramers<2,RelDFFull>>
        compute_full(std::array<std::array<std::shared_ptr<const Matrix>,4>,2> rocoeff, std::array<std::array<std::shared_ptr<const Matrix>,4>,2> iocoeff,
                     std::array<std::list<std::shared_ptr<RelDFHalf>>,2> half, const bool appj, const bool appjj = false);
};


class RelJop : public RelMOFile {
  protected:
    std::shared_ptr<Kramers<2,ZMatrix>> compute_mo1e(std::shared_ptr<const Kramers<2,ZMatrix>> coeff) override;
    std::shared_ptr<Kramers<4,ZMatrix>> compute_mo2e(std::shared_ptr<const Kramers<2,ZMatrix>> coeff) override;

  public:
    RelJop(const std::shared_ptr<const Geometry> geom, const int nstart, const int nfence, std::shared_ptr<const RelCoeff_Block> coeff, const int charge,
      const bool gaunt, const bool breit, const bool restricted = false, const bool tsymm = true)
      : RelMOFile(geom, coeff, charge, gaunt, breit, tsymm) { init(nstart, nfence, restricted); }
};

}

#endif
