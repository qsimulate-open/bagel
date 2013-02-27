//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __SRC_REL_DFHALFCOMPLEX_H
#define __SRC_REL_DFHALFCOMPLEX_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>
#include <src/rel/dfdata.h>
#include <src/rel/alpha.h>
#include <src/rel/breit2index.h>
#include <src/rel/reldfbase.h>
#include <src/df/df.h>

namespace bagel {

class DFData;

class DFHalfComplex : public RelDFBase {
  protected:
    std::array<std::shared_ptr<DFHalfDist>,2> dfhalf_;
    std::array<std::shared_ptr<DFHalfDist>,2> df2_;
    bool split_;

    void set_basis() override;

  public:
    DFHalfComplex(std::shared_ptr<const DFData>, std::vector<std::shared_ptr<ABcases>> bas,
                  std::array<std::shared_ptr<const Matrix>,4>, std::array<std::shared_ptr<const Matrix>,4>);

    DFHalfComplex(std::array<std::shared_ptr<DFHalfDist>,2> data, std::pair<int,int> coord, std::vector<std::shared_ptr<ABcases>> bas);
                  
    std::array<std::shared_ptr<DFHalfDist>, 2> get_data() const { return dfhalf_; }
    std::shared_ptr<DFHalfDist> get_real() const { return dfhalf_[0]; }
    std::shared_ptr<DFHalfDist> get_imag() const { return dfhalf_[1]; }

    bool matches(std::shared_ptr<DFHalfComplex>) const;
    bool alpha_matches(std::shared_ptr<DFHalfComplex>) const;
    bool alpha_matches(std::shared_ptr<Breit2Index>) const;
    std::shared_ptr<DFHalfComplex> multiply_breit(std::shared_ptr<Breit2Index>) const;

    // zaxpy
    void zaxpy(std::complex<double> a, std::shared_ptr<const DFHalfComplex> o);

    // for the zgemm3m-like algorithm
    void set_sum_diff();
    std::shared_ptr<DFHalfDist> sum() const { return df2_[0]; } 
    std::shared_ptr<DFHalfDist> diff() const { return df2_[1]; } 

    std::complex<double> fac() const { assert(basis_.size() == 1); return basis_[0]->fac(); }
    std::list<std::shared_ptr<DFHalfComplex>> split();
    bool split_status() const { return split_; }

};

}

#endif
