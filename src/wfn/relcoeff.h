//
// BAGEL - Parallel electron correlation program.
// Filename: relcoeff.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// 4-component coefficient matrices for use in ZFCI and ZCASSCF
// Basis spinors (rows) are always stored in order L+, L-, S+, S-
// Three formats for MO ordering:
//     Striped format:  A+ A- B+ B- C+ C- D+ D-...
//     Block format:    A+ B+ C+ D+ A- B- C- D-...
//     Kramers format:  Assumes we're coming from Quaternion Diagonalization - on construction, rearranges to block format except with all + before all -
//     For all formats, spaces are stored as occupied, active, virtual, positronic

#ifndef __SRC_WFN_RELCOEFF_H
#define __SRC_WFN_RELCOEFF_H

#include <src/wfn/geometry.h>

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
    int nvirt_nr() const { return nvirt_nr_; }
    int nvirt_rel() const { return nvirt_nr_ + nneg_/2; }

    // spin orbitals (1 column)
    int nneg() const { return nneg_; }
    int npos() const { return 2*(nclosed_ + nact_ + nvirt_nr_); }

    void print_info() const;

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
    RelCoeff_Striped(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false);

    // construct a blank RelCoeff_Striped
    //RelCoeff_Striped(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
    // : RelCoeff(_ndim, _loc, _nclosed, _nact, _nvirt, _nneg) { }

    std::shared_ptr<RelCoeff_Striped> electronic_part() const {
      ZMatrix tmp = slice(0, npos());
      auto out = std::make_shared<RelCoeff_Striped>(tmp, nclosed_, nact_, nvirt_nr_, 0);
      return out;
    }

    std::shared_ptr<RelCoeff_Block> block_format() const;

    // get Kramers-adapted coefficient via quaternion diagonalization
    std::shared_ptr<RelCoeff_Striped> init_kramers_coeff_dirac(std::shared_ptr<const Geometry> geom, std::shared_ptr<const ZMatrix> overlap,
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
    RelCoeff_Block(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg);

    // construct a blank RelCoeff_Block
    //RelCoeff_Block(const int _ndim, const bool _loc, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
    // : RelCoeff(_ndim, _loc, _nclosed, _nact, _nvirt, _nneg) { }

    std::shared_ptr<RelCoeff_Striped> striped_format() const;
};


class RelCoeff_Kramers : public RelCoeff {
  public:
    // standard constructor; copy data directly from _coeff
    RelCoeff_Kramers(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg);
    std::shared_ptr<RelCoeff_Block> block_format() const;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Striped)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Block)

#endif
