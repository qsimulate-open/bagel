//
// BAGEL - Parallel electron correlation program.
// Filename: hso.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_MOLECULE_HSO_H
#define __SRC_MOLECULE_HSO_H

#include <src/math/matrix.h>

namespace bagel {

class HSO {
  protected:
    // Im{aa}, Re{ab}, and Im{ab} components of the spin-orbit hcore
    std::shared_ptr<Matrix> soiaa_, sorab_, soiab_;

  private:
    // serialization
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & soiaa_ & soiab_ & soiab_;
    }

  public:
    HSO() { }
    HSO(const int nbasis) : soiaa_(std::make_shared<Matrix>(nbasis, nbasis)), sorab_(std::make_shared<Matrix>(nbasis, nbasis)),
                            soiab_(std::make_shared<Matrix>(nbasis, nbasis)) {
      soiaa_->zero(); sorab_->zero(); soiab_->zero();
    }
    ~HSO() { }

    void construct_iaa(const int offset1, const int offset0, const int dim1, const int dim0, const double* const data) {
      soiaa_->copy_block(offset1, offset0, dim1, dim0, data);
    }

    void construct_rab(const int offset1, const int offset0, const int dim1, const int dim0, const double* const data) {
      sorab_->copy_block(offset1, offset0, dim1, dim0, data);
    }

    void construct_iab(const int offset1, const int offset0, const int dim1, const int dim0, const double* const data) {
      soiab_->copy_block(offset1, offset0, dim1, dim0, data);
    }

    void allreduce() { // for parallel run
      soiaa_->allreduce();
      sorab_->allreduce();
      soiab_->allreduce();
    }

    void fill_upper_negative() {
      soiaa_->fill_upper_negative();
      sorab_->fill_upper_negative();
      soiab_->fill_upper_negative();
    }

    std::shared_ptr<const Matrix> soiaa() const { return soiaa_;}
    std::shared_ptr<const Matrix> sorab() const { return sorab_;}
    std::shared_ptr<const Matrix> soiab() const { return soiab_;}
};

}

#endif
