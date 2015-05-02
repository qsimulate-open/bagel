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
// Two formats for MO ordering:
//     Striped format:  A+ A- B+ B- C+ C- D+ D-...
//     Block format:    A+ B+ C+ D+ A- B- C- D-...
//     For both formats, spaces are stored as occupied, active, virtual, positronic

#ifndef __SRC_WFN_RELCOEFF_H
#define __SRC_WFN_RELCOEFF_H

#include <src/wfn/geometry.h>

namespace bagel {

class RelCoeff_Striped;
class RelCoeff_Block;

class RelCoeff : public ZMatrix {
  protected:
    RelCoeff(const ZMatrix& _coeff, const int _nclo, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false);

    int nbasis_;
    int nclosed_;
    int nact_;
    int nvirt_nr_;
    int nneg_;

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

    using Matrix_base<std::complex<double>>::copy_block;
};

class RelCoeff_Striped : public RelCoeff {
  protected:

  public:
    RelCoeff_Striped(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false)
   : RelCoeff(_coeff, _nclosed, _nact, _nvirt, _nneg, move_neg) { }

    std::shared_ptr<RelCoeff_Block> block_format() const;
};


class RelCoeff_Block : public RelCoeff {
  protected:

  public:
    RelCoeff_Block(const ZMatrix& _coeff, const int _nclosed, const int _nact, const int _nvirt, const int _nneg)
   : RelCoeff(_coeff, _nclosed, _nact, _nvirt, _nneg) { }

    std::shared_ptr<RelCoeff_Striped> striped_format() const;
};

}

#if 0
#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Striped)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Block)
#endif

#endif
