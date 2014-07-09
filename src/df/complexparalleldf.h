//
// BAGEL - Parallel electron correlation program.
// Filename: complexparalleldf.h
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


#ifndef __SRC_DF_COMPLEXPARALLELDF_H
#define __SRC_DF_COMPLEXPARALLELDF_H

#include <src/math/zmatrix.h>
#include <src/df/paralleldf.h>

namespace bagel {

class ComplexParallelDF : public std::enable_shared_from_this<ComplexParallelDF> {
  protected:

    std::array<std::shared_ptr<ParallelDF>,2> dfdata_;
    const int nblock_;

    // naux runs fastest, nindex2 runs slowest
    const size_t naux_;
    const size_t nindex1_;
    const size_t nindex2_;

    std::shared_ptr<const ComplexParallelDF> df_;
    // data2_ is usually empty (except for the original DFDist)
    // AO two-index integrals ^ -1/2
    std::shared_ptr<Matrix> data2_;

    bool comp_;
    bool swap_;
    bool serial_;

  public:
    ComplexParallelDF(const size_t, const size_t, const size_t, const int, std::shared_ptr<const ComplexParallelDF> = nullptr, std::shared_ptr<Matrix> = nullptr);

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }
    size_t size() const { return naux_*nindex1_*nindex2_; }
    int nblock() const { return nblock_; }

    bool swapped() const { return swap_; }
    std::shared_ptr<const ComplexParallelDF> swap() const;

    const std::shared_ptr<const ComplexParallelDF> df() const { return df_; }
    const std::array<std::shared_ptr<ParallelDF>,2>& dfdata() const { return dfdata_; }
    std::shared_ptr<const ParallelDF> get_real() const { return dfdata_[0]; }
    std::shared_ptr<const ParallelDF> get_imag() const { return dfdata_[1]; }
    std::shared_ptr<const DFBlock> real_block(int i) const { return dfdata_[0]->block_[i]; }
    std::shared_ptr<const DFBlock> imag_block(int i) const { return dfdata_[1]->block_[i]; }

    std::shared_ptr<const StaticDist> adist_now() const { return dfdata_[0]->block_[0]->adist_now(); }

    std::shared_ptr<ZMatrix> form_2index(std::shared_ptr<const ComplexParallelDF> o, const double a, const bool swap = false) const;
    //std::shared_ptr<ZMatrix> form_4index(std::shared_ptr<const ComplexParallelDF> o, const double a, const bool swap = false) const;
    //std::shared_ptr<ZMatrix> form_aux_2index(std::shared_ptr<const ComplexParallelDF> o, const double a) const;

    void ax_plus_y(const double a, const std::shared_ptr<const ComplexParallelDF> o);
    void scale(const double a);
    void symmetrize();

    std::shared_ptr<const Matrix> data2_real() const { return data2_; }
    std::shared_ptr<const ZMatrix> data2() const { return std::make_shared<ZMatrix>(*data2_, 1.0); }

    // compute a J operator, given density matrices in AO basis
    std::shared_ptr<ZMatrix> compute_Jop(const std::shared_ptr<const ZMatrix> den) const;
    std::shared_ptr<ZMatrix> compute_Jop(const std::shared_ptr<const ComplexParallelDF> o, const std::shared_ptr<const ZMatrix> den, const bool onlyonce = false) const;
    std::shared_ptr<ZMatrix> compute_Jop_from_cd(std::shared_ptr<const ZVectorB> cd) const;
    std::shared_ptr<ZVectorB> compute_cd(const std::shared_ptr<const ZMatrix> den, std::shared_ptr<const Matrix> dat2 = nullptr, const bool onlyonce = false) const;
    //std::shared_ptr<ZMatrix> compute_cd(const std::shared_ptr<const Matrix> den, std::shared_ptr<const Matrix> dat2 = nullptr, const bool onlyonce = false) const;

    void average_3index() {
      Timer time;
      if (!serial_)
        for (auto& i : dfdata_[0]->block_)
          i->average();
        for (auto& i : dfdata_[1]->block_)
          i->average();
      time.tick_print("3-index ints post");
    }

    void shell_boundary_3index() {
      throw std::logic_error("Using a function that has not been verified to work with complex integrals - use caution");
      if (!serial_)
        for (auto& i : dfdata_[0]->block_)
          i->shell_boundary();
        for (auto& i : dfdata_[1]->block_)
          i->shell_boundary();
    }
};

}

#endif
