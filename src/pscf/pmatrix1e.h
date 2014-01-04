//
// BAGEL - Parallel electron correlation program.
// Filename: pmatrix1e.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __src_util_pmatrix1e_h
#define __src_util_pmatrix1e_h

#include <vector>
#include <complex>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <src/pscf/pgeometry.h>
#include <src/util/f77.h>
#include <src/pscf/pdata.h>
#include <memory>

namespace bagel {

class PMatrix1e {
  protected:
    std::shared_ptr<PGeometry> geom_;
    const int nbasis_;
    // # of row
    int ndim_;
    // # of column
    int mdim_;
    const int blocksize_;
    const int totalsize_;

    virtual void init();
    virtual void computebatch(const std::vector<std::shared_ptr<const Shell>>&,
        const int, const int, const int, const int) { assert(false); };

    std::shared_ptr<PData> data_;

  public:
    PMatrix1e(const std::shared_ptr<PGeometry>);
    PMatrix1e(const std::shared_ptr<PGeometry>, const int ldn, const int ldm);
    // Constructing PMatrix1e while increasing the leading dimension, or ndim_.
    PMatrix1e(const std::shared_ptr<PMatrix1e> source, const int ldn);
    // Constructing PMatrix1e while reducing the number of columns
    PMatrix1e(const std::shared_ptr<PMatrix1e> source, const std::pair<int, int> mcut);
    // Constructing Pmatrix merging two matrices
    PMatrix1e(const std::shared_ptr<PMatrix1e>, const std::shared_ptr<PMatrix1e>);
    ~PMatrix1e();

    PMatrix1e operator*(const PMatrix1e&) const;
    PMatrix1e operator%(const PMatrix1e&) const; // caution
    PMatrix1e operator+(const PMatrix1e&) const;
    PMatrix1e operator-(const PMatrix1e&) const;
    PMatrix1e& operator=(const PMatrix1e&);
    PMatrix1e& operator+=(const PMatrix1e&);

    PMatrix1e ft() const;
    PMatrix1e bft() const;

    void set_geom(std::shared_ptr<PGeometry> a) {geom_ = a; };

    void hermite();
    void real();
    void conj();
    void conj_transpose();
    void scale(const std::complex<double>);
    const std::complex<double>* bp(const int k) const { return data_->pointer((k + geom_->K()) * blocksize_); };
    std::complex<double>* bpw(const int k) { return data_->pointer((k + geom_->K()) * blocksize_); };

    int nbasis() const { return nbasis_; };
    int K() const { return geom_->K(); };
    int L() const { return geom_->L(); };
    int S() const { return geom_->S(); };
    double A() const { return geom_->A(); };
    std::shared_ptr<PData> data() const { return data_; };
    int ndim() const { return ndim_; };
    int mdim() const { return mdim_; };
    int blocksize() const { return blocksize_; };
    int totalsize() const { return totalsize_; };

    const std::shared_ptr<PGeometry> geom() const { return geom_; };
    void diagonalize(double*);
    std::shared_ptr<PMatrix1e> inverse() const;
    std::tuple<std::shared_ptr<PMatrix1e>, std::shared_ptr<PMatrix1e>> svd();

    void print() const;
    void rprint(const int precision=6) const;

    void ax_plus_y(const std::complex<double>, const PMatrix1e&);
    void ax_plus_y(const std::complex<double>, const std::shared_ptr<PMatrix1e>);
    const std::complex<double> dot_product(const PMatrix1e&) const;
    const std::complex<double> dot_product(const std::shared_ptr<PMatrix1e>) const;

    double rms() const;
    double trace() const;

    std::pair<std::shared_ptr<PMatrix1e>, std::shared_ptr<PMatrix1e>>
      split(const int nrow1, const int nrow2) const;
    std::shared_ptr<PMatrix1e> merge(const std::shared_ptr<PMatrix1e>) const;

    // AO-to-MO integral transformation...
    std::shared_ptr<PMatrix1e> mo_transform(std::shared_ptr<PMatrix1e>,
                             const int istart, const int ifence,
                             const int jstart, const int jfence) const;

};

}

#endif
