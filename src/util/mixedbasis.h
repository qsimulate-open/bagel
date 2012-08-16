//
// Newint - Parallel electron correlation program.
// Filename: mixed_matrix1e.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_UTIL_MIXED_MATRIX1E_H
#define __SRC_UTIL_MIXED_MATRIX1E_H

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cassert>
#include <string>
#include <src/scf/geometry.h>

template<typename T>
class MixedBasis {
  protected:
    const size_t nbasis0_;
    const size_t nbasis1_;

    std::unique_ptr<double[]> data_;

    void computebatch(const std::vector<std::shared_ptr<const Shell> >& input, const int offsetb0, const int offsetb1, const int ndim) {
      // input = [b1, b0]
      assert(input.size() == 2);
      const int dimb1 = input[0]->nbasis();
      const int dimb0 = input[1]->nbasis();
      T batch(input);
      batch.compute();
      const double* odata = batch.data();

      int cnt = 0;
      for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
        for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
          data_[i*ndim + j] = odata[cnt];
        }
      }
    };

  public:
    MixedBasis(const std::shared_ptr<const Geometry> g0, const std::shared_ptr<const Geometry> g1)
     : nbasis0_(g0->nbasis()), nbasis1_(g1->nbasis()), data_(new double[g0->nbasis()*g1->nbasis()]) {
      const std::vector<std::shared_ptr<const Atom> > atoms0 = g0->atoms();
      const std::vector<std::shared_ptr<const Atom> > atoms1 = g1->atoms();
      const std::vector<std::vector<int> > offsets0 = g0->offsets();
      const std::vector<std::vector<int> > offsets1 = g1->offsets();

      // only lower half will be stored
      for (int iatom0 = 0; iatom0 != g0->natom(); ++iatom0) {
        const std::shared_ptr<const Atom> catom0 = atoms0[iatom0];
        const int numshell0 = catom0->shells().size();
        const std::vector<int> coffset0 = offsets0[iatom0];
        const std::vector<std::shared_ptr<const Shell> > shell0 = catom0->shells();

        for (int iatom1 = 0; iatom1 != g1->natom(); ++iatom1) {
          const std::shared_ptr<const Atom> catom1 = atoms1[iatom1];
          const int numshell1 = catom1->shells().size();
          const std::vector<int> coffset1 = offsets1[iatom1];
          const std::vector<std::shared_ptr<const Shell> > shell1 = catom1->shells();

          for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
            const int offset0 = coffset0[ibatch0];
            std::shared_ptr<const Shell> b0 = shell0[ibatch0];
            for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
              const int offset1 = coffset1[ibatch1];
              std::shared_ptr<const Shell> b1 = shell1[ibatch1];
              std::vector<std::shared_ptr<const Shell> > input = {{b1, b0}};

              computebatch(input, offset0, offset1, nbasis1_);

            }
          }
        }
      }

    };
    ~MixedBasis() {};

    const double* data() const { return data_.get(); };
    double* data() { return data_.get(); };

    void print(const std::string in = "", const size_t size = 10) const {
      std::cout << "++++ " << in << " ++++" << std::endl;
      for (int i = 0; i != std::min(size,nbasis1_); ++i) {
        for (int j = 0; j != std::min(size,nbasis0_); ++j) {
          std::cout << std::fixed << std::setw(12) << std::setprecision(9) << data_[j*nbasis1_+i]  << " ";
        }
        std::cout << std::endl;
      }
    };

};

#endif
