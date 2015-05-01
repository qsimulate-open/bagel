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
    int nbasis_;
    int nocc_;
    int nact_;
    int nvirt_;
    int nneg_;

  public:
    RelCoeff(const ZMatrix& _coeff, const int _nocc, const int _nact, const int _nvirt, const int _nneg, const bool move_neg = false);

    int nbasis_nr() const { return nbasis_; }
    int nbasis_rel() const { return 4*nbasis_; }

    // spatial orbitals (2 columns)
    int nocc() const { return nocc_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }

    // spin orbitals (1 column)
    int nneg() const { return nneg_; }
    int npos() const { return 2*(nocc_ + nact_ + nvirt_); }

    using Matrix_base<std::complex<double>>::copy_block;
};

#if 0
class RelCoeff_Striped : public RelCoeff {
  protected:

  public:
    RelCoeff_Striped() { }
    //std::shared_ptr<RelCoeff_Block> block_format();
};


class RelCoeff_Block : public RelCoeff {
  protected:

  public:
    RelCoeff_Block() { }
    //std::shared_ptr<RelCoeff_Striped> striped_format();
};
#endif

}

#if 0
#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Striped)
BOOST_CLASS_EXPORT_KEY(bagel::RelCoeff_Block)
#endif

#endif
