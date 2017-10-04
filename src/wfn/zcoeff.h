//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcoeff.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// Complex coefficient matrices for use in ZFCI and ZCASSCF
// Basis spinors (rows) are always stored in order L+, L-, S+, S-
// Three formats for MO ordering:
//     Striped format:  A+ A- B+ B- C+ C- D+ D-...
//     Block format:    A+ B+ C+ D+ A- B- C- D-... divided as closed, act, {virtual + positronic}
//     Kramers format:  Similar to Block format, but with all + before all -

#ifndef __SRC_WFN_ZCOEFF_H
#define __SRC_WFN_ZCOEFF_H

#include <set>
#include <src/wfn/geometry.h>
#include <src/util/kramers.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class ZCoeff_Striped;
class ZCoeff_Block;

class ZCoeff_base : public ZMatrix {
  protected:
    ZCoeff_base() { }
    ZCoeff_base(const int ndim, const bool loc, const int nclosed, const int nact, const int nvirt, const int nneg);

    int nbasis_;
    int nclosed_;
    int nact_;
    int nvirt_nr_;
    int nneg_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZMatrix>(*this) & nbasis_ & nclosed_ & nact_ & nvirt_nr_ & nneg_;
    }

  public:
    int nbasis_nr() const { return nbasis_; }
    int nbasis_rel() const { return 4*nbasis_; }

    // spatial orbitals (2 columns)
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nocc() const { return nclosed_ + nact_; }
    int nvirt_nr() const { return nvirt_nr_; }
    int nvirt_rel() const { return nvirt_nr_ + nneg_/2; }

    // spin orbitals (1 column)
    int nneg() const { return nneg_; }
    int npos() const { return 2*(nclosed_ + nact_ + nvirt_nr_); }

    void print_info() const;

    // static function used to make the order of eigenvalues & eigenvectors of ZMatrix::diagonalize() match that given by QuatMatrix::diagonalize()
    static void rearrange_eig(VectorB& eig, std::shared_ptr<ZMatrix> coeff, const bool includes_neg = true);

    using Matrix_base<std::complex<double>>::copy_block;
};


class ZCoeff_Striped : public ZCoeff_base {
  private:
    ZCoeff_Striped() { }
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & boost::serialization::base_object<ZCoeff_base>(*this); }

  public:
    // standard constructor; copy data directly from _coeff
    ZCoeff_Striped(const ZMatView& _coeff, const int nclosed, const int nact, const int nvirt, const int nneg = 0, const bool move_neg = false);

    // construct an empty ZCoeff_Striped
    ZCoeff_Striped(const int ndim, const bool loc, const int nclosed, const int nact, const int nvirt, const int nneg = 0)
     : ZCoeff_base(ndim, loc, nclosed, nact, nvirt, nneg) { }

    // copy construct
    ZCoeff_Striped(const ZCoeff_Striped& o) : ZCoeff_Striped(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<ZCoeff_Striped> copy() const { return std::make_shared<ZCoeff_Striped>(*this); }
    std::shared_ptr<ZCoeff_Striped> clone() const { return std::make_shared<ZCoeff_Striped>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<ZCoeff_Striped> electronic_part() const {
      return std::make_shared<ZCoeff_Striped>(slice(0, npos()), nclosed_, nact_, nvirt_nr_, 0);
    }

    std::shared_ptr<ZCoeff_Block> block_format(int nclosed = -1, int nact = -1, int nvirt = -1, int nneg = -1) const;
    std::shared_ptr<Kramers<1,ZMatrix>> kramers_active() const;

    // get Kramers-adapted coefficient via quaternion diagonalization
    std::shared_ptr<const ZCoeff_Striped> init_kramers_coeff(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> overlap,
                                                             std::shared_ptr<const ZMatrix> hcore, const int nele, const bool gaunt, const bool breit) const;

    // rearrange coefficient to {c,a,v} by selecting active columns from input coefficient
    std::shared_ptr<const ZCoeff_Striped> set_active(std::set<int> active_indices, const int nele) const;
};


class ZCoeff_Block : public ZCoeff_base {
  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZCoeff_base>(*this);
    }

  private:
    ZCoeff_Block() { }

  public:
    // standard constructor; copy data directly from _coeff
    ZCoeff_Block(const ZMatView& _coeff, const int nclosed, const int nact, const int nvirt, const int nneg = 0);

    // construct an empty ZCoeff_Block
    ZCoeff_Block(const int ndim, const bool loc, const int nclosed, const int nact, const int nvirt, const int nneg = 0)
     : ZCoeff_base(ndim, loc, nclosed, nact, nvirt, nneg) { }

    // copy construct
    ZCoeff_Block(const ZCoeff_Block& o) : ZCoeff_Block(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<ZCoeff_Block> copy() const { return std::make_shared<ZCoeff_Block>(*this); }
    std::shared_ptr<ZCoeff_Block> clone() const { return std::make_shared<ZCoeff_Block>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<ZCoeff_Block> closed_part() const;
    std::shared_ptr<ZCoeff_Block> active_part() const;
    std::shared_ptr<ZCoeff_Block> electronic_part() const;

    std::shared_ptr<Kramers<1,ZMatrix>> kramers_active() const;
    std::shared_ptr<ZCoeff_Striped> striped_format() const;
};


class ZCoeff_Kramers : public ZCoeff_base {
  private:
    ZCoeff_Kramers() { }
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & boost::serialization::base_object<ZCoeff_base>(*this); }

  public:
    // standard constructor; copy data from _coeff directly
    ZCoeff_Kramers(const ZMatView& _coeff, const int nclosed, const int nact, const int nvirt, const int nneg, const bool move_neg = false);

    // construct an empty ZCoeff_Kramers
    ZCoeff_Kramers(const int ndim, const bool loc, const int nclosed, const int nact, const int nvirt, const int nneg)
     : ZCoeff_base(ndim, loc, nclosed, nact, nvirt, nneg) { }

    // copy construct
    ZCoeff_Kramers(const ZCoeff_Kramers& o) : ZCoeff_Kramers(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<ZCoeff_Kramers> copy() const { return std::make_shared<ZCoeff_Kramers>(*this); }
    std::shared_ptr<ZCoeff_Kramers> clone() const { return std::make_shared<ZCoeff_Kramers>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<ZCoeff_Block> block_format() const;
    std::shared_ptr<ZCoeff_Striped> striped_format() const;

    // For converting from { L+ L- S+ S- } to { L+ S+ L- S- }
    std::shared_ptr<ZCoeff_Kramers> swap_central() const;

    // For converting from { S+ L+ S- L- } to { L+ S+ L- S- } such as after quaternion diagonalization
    std::shared_ptr<ZCoeff_Kramers> move_positronic() const;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff_base)
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff_Striped)
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff_Block)
BOOST_CLASS_EXPORT_KEY(bagel::ZCoeff_Kramers)

#endif
