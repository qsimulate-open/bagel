//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfock.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_REL_DFOCK_H
#define __SRC_REL_DFOCK_H

#include <src/df/reldf.h>

namespace bagel {

class DFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    const bool gaunt_;
    const bool breit_;

    void two_electron_part(const ZMatView coeff, const double scale_ex, const double scale_coulomb);


    void add_Jop_block(std::shared_ptr<const RelDF>, std::list<std::shared_ptr<const RelCDMatrix>>, const double scale);
    void add_Exop_block(std::shared_ptr<const RelDFHalf>, std::shared_ptr<const RelDFHalf>, const double scale, const bool diag = false);
    void driver(std::shared_ptr<const ZMatrix> coeff, bool gaunt, bool breit, const double scale_exchange, const double scale_coulomb);

    // when gradient is requested, we store half-transformed integrals
    // TODO want to avoid "mutable" but this lets us discard integrals later to free up memory
    mutable bool store_half_;
    mutable bool store_half_gaunt_;
    mutable std::list<std::shared_ptr<RelDFHalf>> half_coulomb_;
    mutable std::list<std::shared_ptr<RelDFHalf>> half_gaunt_;
    mutable std::list<std::shared_ptr<RelDFHalf>> half_breit_;

    // if true, do not use bra-ket symmetry in the exchange build (only useful for breit when accurate orbitals are needed).
    bool robust_;

  public:
    DFock(std::shared_ptr<const Geometry> a,  std::shared_ptr<const ZMatrix> hc, const ZMatView coeff, const bool gaunt, const bool breit,
          const bool store_half, const bool robust = false, const double scale_exch = 1.0, const double scale_coulomb = 1.0, const bool store_half_gaunt = false);
    // same as above
    DFock(std::shared_ptr<const Geometry> a, std::shared_ptr<const ZMatrix> hc, std::shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit,
          const bool store_half, const bool robust = false, const double scale_exch = 1.0, const double scale_coulomb = 1.0, const bool store_half_gaunt = false)
     : DFock(a, hc, *coeff, gaunt, breit, store_half, robust, scale_exch, scale_coulomb, store_half_gaunt) {
    }
    // DFock from half-transformed integrals
    DFock(std::shared_ptr<const Geometry> a, std::shared_ptr<const ZMatrix> hc, std::shared_ptr<const ZMatrix> coeff, std::shared_ptr<const ZMatrix> tcoeff,
          std::list<std::shared_ptr<const RelDFHalf>> int1c, std::list<std::shared_ptr<const RelDFHalf>> int2c,
          const double scale_exch = 1.0, const double scale_coulomb = 1.0);

    // Utility functions. They are static so that it could be used from gradient codes

    // T needs to have "zaxpy" and "matches" functions.
    template<class T> static void factorize(T& m) {
      for (auto i = m.begin(); i != m.end(); ++i)
        for (auto j = i; j != m.end(); ) {
          if (i != j && (*i)->matches(*j)) {
            (*i)->ax_plus_y((*j)->fac() / (*i)->fac(), *j);
            j = m.erase(j);
          } else
            ++j;
        }
    }
    static std::list<std::shared_ptr<RelDF>> make_dfdists(std::vector<std::shared_ptr<const DFDist>>, bool);
    static std::list<std::shared_ptr<RelDFHalf>> make_half_complex(std::list<std::shared_ptr<RelDF>>, std::shared_ptr<const ZMatrix>);

    std::list<std::shared_ptr<RelDFHalf>> half_coulomb() const { assert(store_half_); return half_coulomb_; }
    std::list<std::shared_ptr<RelDFHalf>> half_gaunt() const { assert(store_half_gaunt_); return half_gaunt_; }
    std::list<std::shared_ptr<RelDFHalf>> half_breit() const { assert(store_half_gaunt_); return half_breit_; }

    void discard_half() const {
      store_half_ = false;
      store_half_gaunt_ = false;
      half_coulomb_.clear();
      half_gaunt_.clear();
      half_breit_.clear();
    }

    void build_j(std::list<std::shared_ptr<RelDFHalf>> half1, std::list<std::shared_ptr<RelDFHalf>> half2, std::shared_ptr<const ZMatrix> coeff,
                 const bool gaunt, const bool breit, const double scale_coulomb = 1.0, const int number_of_j = 1);
    void build_j(std::list<std::shared_ptr<const RelDFHalf>> half1, std::list<std::shared_ptr<const RelDFHalf>> half2, std::shared_ptr<const ZMatrix> coeff,
                 const bool gaunt, const bool breit, const double scale_coulomb = 1.0, const int number_of_j = 1);
    void build_k(std::list<std::shared_ptr<RelDFHalf>> half1, std::list<std::shared_ptr<RelDFHalf>> half2, std::shared_ptr<const ZMatrix> coeff,
                 const bool gaunt, const bool breit, const double scale_exch = 1.0);
    void build_k(std::list<std::shared_ptr<const RelDFHalf>> half1, std::list<std::shared_ptr<const RelDFHalf>> half2, std::shared_ptr<const ZMatrix> coeff,
                 const bool gaunt, const bool breit, const double scale_exch = 1.0);
};

}

#endif
