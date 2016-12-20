//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: paralleldf.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_DF_PARALLELDF_H
#define __SRC_DF_PARALLELDF_H

#include <src/df/dfinttask_old.h>
#include <src/df/dfinttask.h>
#include <src/df/dfblock.h>

namespace bagel {

class ParallelDF : public std::enable_shared_from_this<ParallelDF> {
  protected:
    // blocks that this process has
    std::vector<std::shared_ptr<DFBlock>> block_;

    // naux runs fastest, nindex2 runs slowest
    const size_t naux_;
    const size_t nindex1_;
    const size_t nindex2_;

    std::shared_ptr<const ParallelDF> df_;
    // data2_ is usually empty (except for the original DFDist)
    // AO two-index integrals ^ -1/2
    std::shared_ptr<Matrix> data2_;

    bool serial_;

  public:
    ParallelDF(const size_t, const size_t, const size_t, std::shared_ptr<const ParallelDF> = nullptr, std::shared_ptr<Matrix> = nullptr, const bool serial = false);
    virtual ~ParallelDF() { }

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }
    size_t size() const { return naux_*nindex1_*nindex2_; }

    bool serial() const { return serial_; }

    std::vector<std::shared_ptr<DFBlock>>& block() { return block_; }
    const std::vector<std::shared_ptr<DFBlock>>& block() const { return block_; }
    std::shared_ptr<DFBlock> block(const size_t i) { return block_[i]; }
    std::shared_ptr<const DFBlock> block(const size_t i) const { return block_[i]; }

    std::shared_ptr<const StaticDist> adist_now() const { return block_[0]->adist_now(); }

    void add_block(std::shared_ptr<DFBlock> o);

    std::shared_ptr<Matrix> form_2index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::shared_ptr<Matrix> form_4index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::shared_ptr<Matrix> form_aux_2index(std::shared_ptr<const ParallelDF> o, const double a) const;

    void add_direct_product(std::shared_ptr<const VectorB> a, std::shared_ptr<const Matrix> b, const double fac)
       { add_direct_product(std::vector<std::shared_ptr<const VectorB>>{a}, std::vector<std::shared_ptr<const Matrix>>{b}, fac); }
    void add_direct_product(std::vector<std::shared_ptr<const VectorB>> a, std::vector<std::shared_ptr<const Matrix>> b, const double fac);

    void ax_plus_y(const double a, const std::shared_ptr<const ParallelDF> o);
    void scale(const double a);
    void symmetrize();

    std::shared_ptr<btas::Tensor3<double>> get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const;

    const std::shared_ptr<const ParallelDF> df() const { return df_; }
    std::shared_ptr<const Matrix> data2() const { return data2_; }

    // compute a J operator, given density matrices in AO basis
    std::shared_ptr<Matrix> compute_Jop(const std::shared_ptr<const Matrix> den) const;
    std::shared_ptr<Matrix> compute_Jop(const std::shared_ptr<const ParallelDF> o, const std::shared_ptr<const Matrix> den, const bool onlyonce = false) const;
    std::shared_ptr<Matrix> compute_Jop_from_cd(std::shared_ptr<const VectorB> cd) const;
    std::shared_ptr<VectorB> compute_cd(const std::shared_ptr<const Matrix> den, std::shared_ptr<const Matrix> dat2 = nullptr, const int number_of_j = 2) const;

    void average_3index() {
      Timer time;
      if (!serial_)
        for (auto& i : block_)
          i->average();
      time.tick_print("3-index ints post");
    }

    void shell_boundary_3index() {
      if (!serial_)
        for (auto& i : block_)
          i->shell_boundary();
    }
};

}

#endif
