//
// BAGEL - Parallel electron correlation program.
// Filename: moint.h
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

// A base class for electorn Correlation methods
// Certain kinds of MO integrals are formed.
//   - aaii (assumes DF - TODO half transformed DF vector might be available..)
//

#ifndef __SRC_SMITH_MOINT_H
#define __SRC_SMITH_MOINT_H

#include <stddef.h>
#include <memory>
#include <stdexcept>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/algo.h>

namespace bagel {
namespace SMITH {

// the template parameter T specifies the storage type

class K2ext {
  protected:
    std::shared_ptr<const SMITH_Info> ref_;
    std::shared_ptr<const Coeff> coeff_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor> data_;

    // some handwritten drivers
    std::map<size_t, std::shared_ptr<DFFullDist>> generate_list();

    void form_4index(const std::map<size_t, std::shared_ptr<DFFullDist>>& dflist);

  public:
    K2ext(std::shared_ptr<const SMITH_Info> r, std::shared_ptr<const Coeff> c, std::vector<IndexRange> b);
    ~K2ext() {}

    std::shared_ptr<Tensor> data() { return data_; }
    std::shared_ptr<Tensor> tensor() { return data_; }

};


class MOFock {
  protected:
    std::shared_ptr<const SMITH_Info> ref_;
    std::shared_ptr<Coeff> coeff_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor> data_;
    std::shared_ptr<Tensor> hcore_;

  public:
    MOFock(std::shared_ptr<const SMITH_Info> r, std::vector<IndexRange> b);
    ~MOFock() {}

    std::shared_ptr<Tensor> data() { return data_; }
    std::shared_ptr<Tensor> tensor() { return data_; }
    std::shared_ptr<Tensor> hcore() { return hcore_; }
    std::shared_ptr<const Coeff> coeff() const { return coeff_; }
};


class Ci {
  protected:
    std::shared_ptr<const SMITH_Info> ref_;
    std::vector<IndexRange> blocks_;
    std::size_t ci_size_;
    std::shared_ptr<Tensor>  rdm0deriv_;


  public:
    Ci(std::shared_ptr<const SMITH_Info> r, std::vector<IndexRange> b, std::shared_ptr<const Civec> c);
    ~Ci() {}

    std::shared_ptr<Tensor> tensor() const { return rdm0deriv_; }
};

}
}

#endif

