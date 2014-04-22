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

#include <numeric>
#include <src/math/zmatrix.h>
#include <src/util/timer.h>
#include <src/parallel/staticdist.h>

namespace bagel {

/*
    DFBlock_London is a slice of 3-index DF integrals. Distributed by the first index
*/

class DFBlock_London {
  protected:
    // aux_ runs fastest, b2_ runs slowest
    std::unique_ptr<std::complex<double>[]> data_;

    // distribution information
    const std::shared_ptr<const StaticDist> adist_shell_;
    const std::shared_ptr<const StaticDist> adist_;
    // if true, asize is evenly distributed. If false, asize is at the shell boundary
    bool averaged_;

    // dimensions of the block
    size_t asize_;
    size_t b1size_;
    size_t b2size_;

    // a set of offsets of this block in the entire DF integrals
    size_t astart_;
    size_t b1start_;
    size_t b2start_;

  public:
    // construction of a block from AO integrals
    DFBlock_London(std::shared_ptr<const StaticDist> adist_shell, std::shared_ptr<const StaticDist> adist,
            const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s, const bool averaged = false);

    DFBlock_London(const DFBlock_London& o);

    // average the asize between MPI processes (block will be described by dist_)
    void average();
    void shell_boundary();

    std::shared_ptr<DFBlock_London> transform_second(std::shared_ptr<const ZMatrix> c, const bool trans = false) const;
    std::shared_ptr<DFBlock_London> transform_third(std::shared_ptr<const ZMatrix> c, const bool trans = false) const;

    std::shared_ptr<DFBlock_London> clone() const;
    std::shared_ptr<DFBlock_London> copy() const;
    void zero() { std::fill_n(data_.get(), size(), 0.0); }

    // dist
    const std::shared_ptr<const StaticDist>& adist_now() const { return averaged_ ? adist_ : adist_shell_; }

    // dimensions of the block
    size_t asize() const { return asize_; }
    size_t b1size() const { return b1size_; }
    size_t b2size() const { return b2size_; }

    size_t size() const { return asize_*b1size_*b2size_; }
    bool averaged() const { return averaged_; }

    // a set of offsets of this block in the entire DF integrals
    size_t astart() const { return astart_; }
    size_t b1start() const { return b1start_; }
    size_t b2start() const { return b2start_; }

    // TODO direct access to data will be disabled once implementation is done
    std::complex<double>* get() { return data_.get(); }
    const std::complex<double>* get() const { return data_.get(); }
    std::complex<double>& operator[](const size_t i) { return data_[i]; }

    // some math functions
    DFBlock_London& operator+=(const DFBlock_London& o);
    DFBlock_London& operator-=(const DFBlock_London& o);
    void ax_plus_y(const double a, const DFBlock_London& o);
    void ax_plus_y(const double a, const std::shared_ptr<const DFBlock_London> o) { ax_plus_y(a, *o); }
    void scale(const double a);

    // add ab^+  to this.
    void add_direct_product(const std::shared_ptr<const ZMatrix> a, const std::shared_ptr<const ZMatrix> b, const double fac);

    // some additional functions
    // symmetrize b1 and b2 (assuming b1size_ == b2size_)
    void symmetrize();

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

    std::unique_ptr<std::complex<double>[]> form_vec(const std::shared_ptr<const ZMatrix> den) const;
    std::shared_ptr<ZMatrix> form_mat(const std::complex<double>* fit) const;

    void contrib_apply_J(const std::shared_ptr<const DFBlock_London> o, const std::shared_ptr<const ZMatrix> mat);

    void copy_block(const std::shared_ptr<const ZMatrix> o, const int jdim, const size_t offset);
    void add_block(const std::shared_ptr<const ZMatrix> o, const int jdim, const size_t offset, const double fac = 1.0);
    // compute (D|ia)(ia|j) and set to the location specified by the offset
    std::shared_ptr<ZMatrix> form_Dj(const std::shared_ptr<const ZMatrix> o, const int jdim) const;

    // CAUTION, ist, jst, and kst are absolute number (NOT relative to astart_, ...). Returns complex<double>[] whose size is i*j*k
    std::shared_ptr<ZMatrix> get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

    // get_block returns (j|cd), while get_block_conj returns (j*|cd)
    std::shared_ptr<ZMatrix> get_block_conj(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

    // use with caution
    std::unique_ptr<std::complex<double>[]> release_data() { return std::move(data_); }
};

}


#endif
