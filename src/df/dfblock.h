//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_DF_DFBLOCK_H
#define __SRC_DF_DFBLOCK_H

#include <src/math/matrix.h>
#include <src/math/matop.h>
#include <src/df/dfblock_base.h>

namespace bagel {

/*
    DFBlock is a slice of 3-index DF integrals. Distributed by the first index
*/

class DFBlock : public DFBlock_base<double> {
  public:
    template<typename... Types>
    DFBlock(Types&&... args) : DFBlock_base<double>(std::forward<Types>(args)...) { }

    DFBlock(const DFBlock& o) : DFBlock_base<double>(o) { }
    DFBlock(DFBlock&& o)      : DFBlock_base<double>(std::move(o)) { }

    DFBlock& operator=(const DFBlock& o) { DFBlock_base<double>::operator=(o); return *this; }
    DFBlock& operator=(DFBlock&& o)      { DFBlock_base<double>::operator=(std::move(o)); return *this; }
    DFBlock& operator+=(const DFBlock& o){ DFBlock_base<double>::operator+=(o); return *this; }
    DFBlock& operator-=(const DFBlock& o){ DFBlock_base<double>::operator-=(o); return *this; }

    std::shared_ptr<DFBlock> clone() const;
    std::shared_ptr<DFBlock> copy() const;

    std::shared_ptr<DFBlock> transform_second(std::shared_ptr<const MatView> c, const bool trans = false) const;
    std::shared_ptr<DFBlock> transform_third(std::shared_ptr<const MatView> c, const bool trans = false) const;
    // TODO will be deprecated
    std::shared_ptr<DFBlock> transform_second(std::shared_ptr<const Matrix> c, const bool trans = false) const;
    std::shared_ptr<DFBlock> transform_third(std::shared_ptr<const Matrix> c, const bool trans = false) const;

    // add ab^+  to this.
    void add_direct_product(const std::shared_ptr<const Matrix> a, const std::shared_ptr<const Matrix> b, const double fac);

    // exchange b1 and b2
    std::shared_ptr<DFBlock> swap() const;

    // 2RDM contractions
    std::shared_ptr<DFBlock> apply_rhf_2RDM(const double scale_exch) const;
    std::shared_ptr<DFBlock> apply_uhf_2RDM(const btas::Tensor2<double>&, const btas::Tensor2<double>&) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm) const;

    // Form 2- and 4-index integrals
    std::shared_ptr<Matrix> form_2index(const std::shared_ptr<const DFBlock> o, const double a) const;
    std::shared_ptr<Matrix> form_4index(const std::shared_ptr<const DFBlock> o, const double a) const;
    // slowest index of o is fixed to n
    std::shared_ptr<Matrix> form_4index_1fixed(const std::shared_ptr<const DFBlock> o, const double a, const size_t n) const;
    std::shared_ptr<Matrix> form_aux_2index(const std::shared_ptr<const DFBlock> o, const double a) const;

    std::unique_ptr<double[]> form_vec(const std::shared_ptr<const Matrix> den) const;
    std::shared_ptr<Matrix> form_mat(const double* fit) const;

    void contrib_apply_J(const std::shared_ptr<const DFBlock> o, const std::shared_ptr<const Matrix> mat);

    // compute (D|ia)(ia|j) and set to the location specified by the offset
    std::shared_ptr<Matrix> form_Dj(const std::shared_ptr<const Matrix> o, const int jdim) const;

    // CAUTION, ist, jst, and kst are absolute number (NOT relative to astart_, ...). Returns double[] whose size is i*j*k
    std::shared_ptr<btas::Tensor3<double>> get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

};

}


#endif
