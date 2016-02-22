//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relcoeff.h
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

// 4-component coefficient matrices for use in ZFCI and ZCASSCF
// Basis spinors (rows) are always stored in order L+, L-, S+, S-
// Three formats for MO ordering:
//     Striped format:  A+ A- B+ B- C+ C- D+ D-...
//     Block format:    A+ B+ C+ D+ A- B- C- D-... divided as closed, act, {virtual + positronic}
//     Kramers format:  Similar to Block format, but with all + before all -

#ifndef __SRC_WFN_RELCOEFF_H
#define __SRC_WFN_RELCOEFF_H

#include <set>
#include <src/wfn/geometry.h>
#include <src/util/kramers.h>

namespace bagel {

class RelCoeff_Striped;
class RelCoeff_Block;

class RelCoeff : public ZMatrix {
  protected:
    RelCoeff(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg);
    RelCoeff() { }

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


class RelCoeff_Striped : public RelCoeff {
  private:
    RelCoeff_Striped() { }
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & boost::serialization::base_object<RelCoeff>(*this); }

  public:
    // standard constructor; copy data directly from _coeff
    RelCoeff_Striped(const ZMatView& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false);

    // construct an empty RelCoeff_Striped
    RelCoeff_Striped(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
     : RelCoeff(_ndim, _loc, _nclosed, _nact, _nvirt, _nneg) { }

    // copy construct
    RelCoeff_Striped(const RelCoeff_Striped& o) : RelCoeff_Striped(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<RelCoeff_Striped> copy() const { return std::make_shared<RelCoeff_Striped>(*this); }
    std::shared_ptr<RelCoeff_Striped> clone() const { return std::make_shared<RelCoeff_Striped>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<RelCoeff_Striped> electronic_part() const {
      return std::make_shared<RelCoeff_Striped>(slice(0, npos()), nclosed_, nact_, nvirt_nr_, 0);
    }

    std::shared_ptr<RelCoeff_Block> block_format(int nclosed = -1, int nact = -1, int nvirt = -1, int nneg = -1) const;
    std::shared_ptr<Kramers<1,ZMatrix>> kramers_active() const;

    // function to generate modified virtual MOs from either a Fock matrix or the one-electron Hamiltonian
    std::shared_ptr<const RelCoeff_Striped> generate_mvo(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> overlap,
              std::shared_ptr<const ZMatrix> hcore, const int ncore, const int nocc_mvo, const bool hcore_mvo, const bool tsymm, const bool gaunt, const bool breit) const;

    // get Kramers-adapted coefficient via quaternion diagonalization
    std::shared_ptr<const RelCoeff_Striped> init_kramers_coeff(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> overlap,
                                                 std::shared_ptr<const ZMatrix> hcore, const int nele, const bool tsymm, const bool gaunt, const bool breit) const;

    // rearrange coefficient to {c,a,v} by selecting active columns from input coefficient
    std::shared_ptr<const RelCoeff_Striped> set_active(std::set<int> active_indices, const int nele, const bool paired) const;
};


class RelCoeff_Block : public RelCoeff {
  private:
    RelCoeff_Block() { }
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & boost::serialization::base_object<RelCoeff>(*this); }

  public:
    // standard constructor; copy data directly from _coeff
    RelCoeff_Block(const ZMatView& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg);

    // construct an empty RelCoeff_Block
    RelCoeff_Block(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
     : RelCoeff(_ndim, _loc, _nclosed, _nact, _nvirt, _nneg) { }

    // copy construct
    RelCoeff_Block(const RelCoeff_Block& o) : RelCoeff_Block(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<RelCoeff_Block> copy() const { return std::make_shared<RelCoeff_Block>(*this); }
    std::shared_ptr<RelCoeff_Block> clone() const { return std::make_shared<RelCoeff_Block>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<RelCoeff_Block> closed_part() const;
    std::shared_ptr<RelCoeff_Block> active_part() const;
    std::shared_ptr<RelCoeff_Block> electronic_part() const;
    std::shared_ptr<RelCoeff_Block> closed_act_positronic() const;
    std::shared_ptr<RelCoeff_Block> update_electronic(std::shared_ptr<const ZMatrix> newcoeff) const;
    std::shared_ptr<RelCoeff_Block> update_closed_act_positronic(std::shared_ptr<const ZMatrix> newcoeff) const;

    std::shared_ptr<RelCoeff_Striped> striped_format() const;
    std::shared_ptr<Kramers<1,ZMatrix>> kramers_active() const;
};


class RelCoeff_Kramers : public RelCoeff {
  private:
    RelCoeff_Kramers() { }
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) { ar & boost::serialization::base_object<RelCoeff>(*this); }

  public:
    // standard constructor; copy data from _coeff directly
    RelCoeff_Kramers(const ZMatView& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false);

    // construct an empty RelCoeff_Kramers
    RelCoeff_Kramers(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
     : RelCoeff(_ndim, _loc, _nclosed, _nact, _nvirt, _nneg) { }

    // copy construct
    RelCoeff_Kramers(const RelCoeff_Kramers& o) : RelCoeff_Kramers(o, o.nclosed_, o.nact_, o.nvirt_nr_, o.nneg_) { }

    std::shared_ptr<RelCoeff_Kramers> copy() const { return std::make_shared<RelCoeff_Kramers>(*this); }
    std::shared_ptr<RelCoeff_Kramers> clone() const { return std::make_shared<RelCoeff_Kramers>(ndim(), localized_, nclosed_, nact_, nvirt_nr_, nneg_); }

    std::shared_ptr<RelCoeff_Block> block_format() const;
    std::shared_ptr<RelCoeff_Striped> striped_format() const;

    // For converting from { L+ L- S+ S- } to { L+ S+ L- S- }
    std::shared_ptr<RelCoeff_Kramers> swap_central() const;

    // For converting from { S+ L+ S- L- } to { L+ S+ L- S- } such as after quaternion diagonalization
    std::shared_ptr<RelCoeff_Kramers> move_positronic() const;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Striped)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Block)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Kramers)

#endif
