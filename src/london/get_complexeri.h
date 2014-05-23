//
// BAGEL - Parallel electron correlation program.
// Filename: get_complexeri.h
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


#ifndef __BAGEL_SRC_LONDON_GET_COMPLEXERI_H
#define __BAGEL_SRC_LONDON_GET_COMPLEXERI_H

#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/wfn/geometry_london.h>
#include <src/math/zmatrix.h>
#include <src/math/matrix.h>
#include <memory>
#include <vector>
#include <complex>

namespace bagel {


std::shared_ptr<ZMatrix> get_ERI (std::shared_ptr<const Geometry_London> cgeom_) {

  const std::vector<std::shared_ptr<const Atom>> atoms = cgeom_->atoms();
  std::vector<std::shared_ptr<const Shell>> basis;
  std::vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const std::vector<std::shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const std::vector<int> tmpoff = cgeom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }

  const int size = basis.size();
  std::cout << "Total number of shells = "  << size << std::endl;
  int nbasis = 0;
  for (int i=0; i!=size; i++) nbasis += basis[i]->nbasis();
  std::cout << "Total number of (contracted) basis functions = "  << nbasis << std::endl;

  const std::shared_ptr<ZMatrix> out = std::make_shared<ZMatrix>(nbasis*nbasis, nbasis*nbasis);
  std::complex<double>* data = out->data();

  for (int i0 = 0; i0 != size; ++i0) {
    const std::shared_ptr<const Shell>  b0 = basis[i0];

    const int b0offset = offset[i0];
    const int b0size = b0->nbasis();
    //std::cout << "i0 = " << i0 << ", b0offset = " << b0offset << ", b0size = " << b0size << std::endl;

    for (int i1 = 0; i1 != size; ++i1) {
      const std::shared_ptr<const Shell>  b1 = basis[i1];

      const int b1offset = offset[i1];
      const int b1size = b1->nbasis();
      //std::cout << "  i1 = " << i1 << ", b1offset = " << b1offset << ", b1size = " << b1size << std::endl;

      for (int i2 = 0; i2 != size; ++i2) {
        const std::shared_ptr<const Shell>  b2 = basis[i2];

        const int b2offset = offset[i2];
        const int b2size = b2->nbasis();
        //std::cout << "    i2 = " << i2 << ", b2offset = " << b2offset << ", b2size = " << b2size << std::endl;

        for (int i3 = 0; i3 != size; ++i3) {
          const std::shared_ptr<const Shell>  b3 = basis[i3];

          const int b3offset = offset[i3];
          const int b3size = b3->nbasis();
          //std::cout << "      i3 = " << i3 << ", b3offset = " << b3offset << ", b3size = " << b3size << std::endl;

          /******************************/
          std::array<std::shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
          ComplexERIBatch eribatch (input, 0.0);
          eribatch.compute();
          const std::complex<double>* eridata = eribatch.data();
          /******************************/

          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int p0 = j0*nbasis;

            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const int p1 = (j1 + p0)*nbasis;

              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int p2 = (j2 + p1)*nbasis;

                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3) {
                  const int p3 = j3 + p2;

                  data[p3] = *eridata;
                  eridata++;

                }
              }
            }
          }
        }
      }
    }
  }
  return out;
}


std::shared_ptr<Matrix> get_ERI_original (std::shared_ptr<const Geometry_London> cgeom_) {

  const std::vector<std::shared_ptr<const Atom>> atoms = cgeom_->atoms();
  std::vector<std::shared_ptr<const Shell>> basis;
  std::vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const std::vector<std::shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const std::vector<int> tmpoff = cgeom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }

  const int size = basis.size();
  std::cout << "Total number of shells = "  << size << std::endl;
  int nbasis = 0;
  for (int i=0; i!=size; i++) nbasis += basis[i]->nbasis();
  std::cout << "Total number of (contracted) basis functions = "  << nbasis << std::endl;

  const std::shared_ptr<Matrix> out = std::make_shared<Matrix>(nbasis*nbasis, nbasis*nbasis);
  double* data = out->data();

  for (int i0 = 0; i0 != size; ++i0) {
    const std::shared_ptr<const Shell>  b0 = basis[i0];

    const int b0offset = offset[i0];
    const int b0size = b0->nbasis();
    //std::cout << "i0 = " << i0 << ", b0offset = " << b0offset << ", b0size = " << b0size << std::endl;

    for (int i1 = 0; i1 != size; ++i1) {
      const std::shared_ptr<const Shell>  b1 = basis[i1];

      const int b1offset = offset[i1];
      const int b1size = b1->nbasis();
      //std::cout << "  i1 = " << i1 << ", b1offset = " << b1offset << ", b1size = " << b1size << std::endl;

      for (int i2 = 0; i2 != size; ++i2) {
        const std::shared_ptr<const Shell>  b2 = basis[i2];

        const int b2offset = offset[i2];
        const int b2size = b2->nbasis();
        //std::cout << "    i2 = " << i2 << ", b2offset = " << b2offset << ", b2size = " << b2size << std::endl;

        for (int i3 = 0; i3 != size; ++i3) {
          const std::shared_ptr<const Shell>  b3 = basis[i3];

          const int b3offset = offset[i3];
          const int b3size = b3->nbasis();
          //std::cout << "      i3 = " << i3 << ", b3offset = " << b3offset << ", b3size = " << b3size << std::endl;

          /******************************/
          std::array<std::shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
          ERIBatch eribatch (input, 0.0);
          eribatch.compute();
          const double* eridata = eribatch.data();
          /******************************/

          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int p0 = j0*nbasis;

            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const int p1 = (j1 + p0)*nbasis;

              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int p2 = (j2 + p1)*nbasis;

                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3) {
                  const int p3 = j3 + p2;

                  data[p3] = *eridata;
                  eridata++;

                }
              }
            }
          }
        }
      }
    }
  }
  return out;
}


}

#endif
