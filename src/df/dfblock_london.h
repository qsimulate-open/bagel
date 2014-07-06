//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock_london.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_DF_DFBLOCK_LONDON_H
#define __SRC_DF_DFBLOCK_LONDON_H

#include <src/math/zmatrix.h>
#include <src/df/dfblock_base.h>
#include <src/math/matop.h>

namespace bagel {

/*
    DFBlock_London is a slice of 3-index DF integrals. Distributed by the first index
*/

class DFBlock_London : public DFBlock_base<std::complex<double>> {
  public:
    template<typename... Types>
    DFBlock_London(Types&&... args) : DFBlock_base<std::complex<double>>(std::forward<Types>(args)...) { }

    DFBlock_London(const DFBlock_London& o) : DFBlock_base<std::complex<double>>(o) { }
    DFBlock_London(DFBlock_London&& o)      : DFBlock_base<std::complex<double>>(std::move(o)) { }

    DFBlock_London& operator=(const DFBlock_London& o) { DFBlock_base<std::complex<double>>::operator=(o); return *this; }
    DFBlock_London& operator=(DFBlock_London&& o)      { DFBlock_base<std::complex<double>>::operator=(std::move(o)); return *this; }
    DFBlock_London& operator+=(const DFBlock_London& o){ DFBlock_base<std::complex<double>>::operator+=(o); return *this; }
    DFBlock_London& operator-=(const DFBlock_London& o){ DFBlock_base<std::complex<double>>::operator-=(o); return *this; }

    std::shared_ptr<DFBlock_London> transform_second(std::shared_ptr<const ZMatView> c, const bool trans = false) const;
    std::shared_ptr<DFBlock_London> transform_third(std::shared_ptr<const ZMatView> c, const bool trans = false) const;
    // TODO will be deprecated
    std::shared_ptr<DFBlock_London> transform_second(std::shared_ptr<const ZMatrix> c, const bool trans = false) const;
    std::shared_ptr<DFBlock_London> transform_third(std::shared_ptr<const ZMatrix> c, const bool trans = false) const;

    std::shared_ptr<DFBlock_London> clone() const;
    std::shared_ptr<DFBlock_London> copy() const;

    // add ab^+  to this.
    void add_direct_product(const std::shared_ptr<const ZMatrix> a, const std::shared_ptr<const ZMatrix> b, const double fac);

    // exchange b1 and b2
    std::shared_ptr<DFBlock_London> swap() const;

/*
    // 2RDM contractions
    std::shared_ptr<DFBlock_London> apply_rhf_2RDM(const double scale_exch) const;
    std::shared_ptr<DFBlock_London> apply_uhf_2RDM(const double*, const double*) const;
    std::shared_ptr<DFBlock_London> apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const;
    std::shared_ptr<DFBlock_London> apply_2RDM(const double* rdm) const;
*/
    // Form 2- and 4-index integrals
    std::shared_ptr<ZMatrix> form_2index(const std::shared_ptr<const DFBlock_London> o, const double a) const;
    std::shared_ptr<ZMatrix> form_4index(const std::shared_ptr<const DFBlock_London> o, const double a) const;
    // slowest index of o is fixed to n
    std::shared_ptr<ZMatrix> form_4index_1fixed(const std::shared_ptr<const DFBlock_London> o, const double a, const size_t n) const;
    std::shared_ptr<ZMatrix> form_aux_2index(const std::shared_ptr<const DFBlock_London> o, const double a) const;

    std::shared_ptr<ZVectorB> form_vec(const std::shared_ptr<const ZMatrix> den) const;
    std::shared_ptr<ZMatrix> form_mat(const std::complex<double>* fit) const;

    void contrib_apply_J(const std::shared_ptr<const DFBlock_London> o, const std::shared_ptr<const ZMatrix> mat);

    // compute (D|ia)(ia|j) and set to the location specified by the offset
    std::shared_ptr<ZMatrix> form_Dj(const std::shared_ptr<const ZMatrix> o, const int jdim) const;

    // CAUTION, ist, jst, and kst are absolute number (NOT relative to astart_, ...). Returns complex<double>[] whose size is i*j*k
    std::shared_ptr<btas::Tensor3<std::complex<double>>> get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

    // get_block returns (j|cd), while get_block_conj returns (j*|cd)
    std::shared_ptr<btas::Tensor3<std::complex<double>>> get_block_conj(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

};

}


#endif
