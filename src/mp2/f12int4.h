//
// Newint - Parallel electron correlation program.
// Filename: f12int4.h
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


// carbon copy of what I wrote in the orz package 
// meant to be standalone

#ifndef __SRC_MP2_F12INT4_H
#define __SRC_MP2_F12INT4_H

#include <src/rysint/rysint.h>
#include <src/rysint/eribatch.h>
#include <src/slater/slaterbatch.h>
#include <src/scf/geometry.h>
#include <src/wfn/reference.h>
#include <src/smith/prim_op.h>
#include <src/mp2/f12mat.h>
#include <tuple>


// Reference implementation of an F12 theory
class F12Ref {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    const int ncore_;
    const double gamma_;
  
  public:
    F12Ref(std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r, const int i, const double gamma)
        : geom_(g), ref_(r), ncore_(i), gamma_(gamma) {};
    ~F12Ref() {};

    void compute();

    std::tuple<std::shared_ptr<Matrix1e>, std::shared_ptr<Matrix1e>, std::shared_ptr<Matrix1e>, int> generate_cabs() const;
};


class File2 {
  protected:
    std::unique_ptr<double[]> data_;
    const size_t dim0_;
    const size_t dim1_;
    const size_t occ_;
  public:
    File2(std::unique_ptr<double[]>& a, const size_t d, const size_t d1, const size_t o) : data_(move(a)), dim0_(d), dim1_(d1), occ_(o) {};
    ~File2() {};

    std::shared_ptr<F12Mat> f12mat(const double* const coeff) const {
      return f12mat(coeff, coeff);
    };
    std::shared_ptr<F12Mat> f12mat(const double* const coeff, const double* const coeff2) const {
      std::unique_ptr<double[]> tmp(new double[occ_*dim1_*occ_*occ_]);
      dgemm_("T", "N", occ_, dim1_*occ_*occ_, dim0_, 1.0, coeff, dim0_, data_.get(), dim0_, 0.0, tmp.get(), occ_);
      std::unique_ptr<double[]> tmp2(new double[occ_*occ_*occ_*occ_]);
      for (size_t i = 0 ; i != occ_*occ_; ++i)
        dgemm_("N", "N", occ_, occ_, dim1_, 1.0, tmp.get()+dim1_*occ_*i, occ_, coeff2, dim1_, 0.0, tmp2.get()+occ_*occ_*i, occ_);
      return std::shared_ptr<F12Mat>(new F12Mat(occ_, std::move(tmp2)));
    };
    std::shared_ptr<F12Ten> f12ten(const double* const coeff, const double* const coeff2, const size_t n, const size_t n2) const {
      std::unique_ptr<double[]> tmp(new double[n*dim1_*occ_*occ_]);
      dgemm_("T", "N", n, dim1_*occ_*occ_, dim0_, 1.0, coeff, dim0_, data_.get(), dim0_, 0.0, tmp.get(), n);
      std::unique_ptr<double[]> tmp2(new double[n*n2*occ_*occ_]);
      for (size_t i = 0 ; i != occ_*occ_; ++i)
        dgemm_("N", "N", n, n2, dim1_, 1.0, tmp.get()+dim1_*n*i, n, coeff2, dim1_, 0.0, tmp2.get()+n*n2*i, n);
      return std::shared_ptr<F12Ten>(new F12Ten(occ_, n, n2, std::move(tmp2))); 
    };
};


class File4 {
  protected:
    std::unique_ptr<double[]> data_;
    const size_t dim_;
    const size_t dim0_;
    const size_t dim1_;
  public:
    File4(std::unique_ptr<double[]>& a, const size_t d, const size_t d0, const size_t d1)
    : data_(move(a)), dim_(d), dim0_(d0), dim1_(d1) {};
    ~File4() {};

    std::shared_ptr<File2> half_transform(const double* const c, const size_t n) const {
      // first transform
      std::unique_ptr<double[]> tmp2(new double[n*dim_*dim0_*dim1_]);
      {
        std::unique_ptr<double[]> tmp(new double[n*dim_*dim0_*dim1_]);
        dgemm_("N", "N", dim0_*dim_*dim1_, n, dim_, 1.0, data_.get(), dim0_*dim_*dim1_, c, dim_, 0.0, tmp.get(), dim0_*dim_*dim1_); 
        // then sort the indices
        SMITH::sort_indices<0,2,1,3,0,1,1,1>(tmp, tmp2, dim0_, dim_, dim1_, n);
      }
      std::unique_ptr<double[]> tmp(new double[n*n*dim0_*dim1_]);
        for (size_t i = 0; i != n; ++i)
          dgemm_("N", "N", dim0_*dim1_, n, dim_, 1.0, tmp2.get()+dim0_*dim1_*dim_*i, dim0_*dim1_, c, dim_, 0.0, tmp.get()+dim0_*dim1_*n*i, dim0_*dim1_); 
      return std::shared_ptr<File2>(new File2(tmp, dim0_, dim1_, n));
    };
};

// unpacked AO integrals
template<typename T> class AOInt {
  protected:
    std::shared_ptr<const Geometry> geom_;
    const size_t dim_;
    const size_t dim0_;
    const size_t dim1_;
    const double gamma_;
    const bool aux0_;
    const bool aux1_;
    std::unique_ptr<double[]> data_;
    std::unique_ptr<double[]> data2_;

    const bool yukawa_;

    void init() {
      // initializing shell info
      std::vector<std::shared_ptr<const Shell> > basis; 
      std::vector<int> offset;
      const std::vector<std::shared_ptr<Atom> > atoms = geom_->atoms(); 
      int cnt = 0;
      for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
        const std::vector<std::shared_ptr<const Shell> > tmp = (*aiter)->shells();
        basis.insert(basis.end(), tmp.begin(), tmp.end());  
        const std::vector<int> tmpoff = geom_->offset(cnt); 
        offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
      }
      // initializing shell info
      std::vector<std::shared_ptr<const Shell> > basis0; 
      std::vector<int> offset0;
      if (aux0_) {
        const std::vector<std::shared_ptr<Atom> > atoms = geom_->aux_atoms(); 
        cnt = 0;
        for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
          const std::vector<std::shared_ptr<const Shell> > tmp = (*aiter)->shells();
          basis0.insert(basis0.end(), tmp.begin(), tmp.end());  
          const std::vector<int> tmpoff = geom_->aux_offset(cnt);
          offset0.insert(offset0.end(), tmpoff.begin(), tmpoff.end());
        }
      } else {
        basis0 = basis;   
        offset0 = offset; 
      }
      // initializing shell info
      std::vector<std::shared_ptr<const Shell> > basis1;
      std::vector<int> offset1;
      if (aux1_) {
        const std::vector<std::shared_ptr<Atom> > atoms = geom_->aux_atoms(); 
        cnt = 0;
        for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
          const std::vector<std::shared_ptr<const Shell> > tmp = (*aiter)->shells();
          basis1.insert(basis1.end(), tmp.begin(), tmp.end());  
          const std::vector<int> tmpoff = geom_->aux_offset(cnt); 
          offset1.insert(offset1.end(), tmpoff.begin(), tmpoff.end());
        }
      } else {
        basis1 = basis;
        offset1 = offset;
      }

      // we will compute redundant integrals as well.
      const size_t size = basis.size();
      for (int i0 = 0; i0 != size; ++i0) {
        const std::shared_ptr<const Shell>  b0 = basis[i0];
        const int b0offset = offset[i0]; 
        const int b0size = b0->nbasis();
        for (int i1 = 0; i1 != basis1.size(); ++i1) {
          const std::shared_ptr<const Shell>  b1 = basis1[i1];
          const int b1offset = offset1[i1]; 
          const int b1size = b1->nbasis();
          for (int i2 = 0; i2 != size; ++i2) {
            const std::shared_ptr<const Shell>  b2 = basis[i2];
            const int b2offset = offset[i2]; 
            const int b2size = b2->nbasis();
            for (int i3 = 0; i3 != basis0.size(); ++i3) {
              const std::shared_ptr<const Shell>  b3 = basis0[i3];
              const int b3offset = offset0[i3]; 
              const int b3size = b3->nbasis();
  
              std::vector<std::shared_ptr<const Shell> > input = {{b3, b2, b1, b0}};
  
              T batch(input, 0.0, gamma_, yukawa_);
              batch.compute();
              const double* eridata = batch.data();
              assert((int)batch.data_size() == b0size * b1size * b2size * b3size);
  
              for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
                for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
                  for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {  
                    std::copy(eridata, eridata+b3size, data_.get()+b3offset+dim0_*(j2+dim_*(j1+dim1_*j0)));
                    eridata += b3size;
                  }
                }
              }
              if (yukawa_) {
                eridata = batch.data2();
                for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
                  for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
                    for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {  
                      std::copy(eridata, eridata+b3size, data2_.get()+b3offset+dim0_*(j2+dim_*(j1+dim1_*j0)));
                      eridata += b3size;
                    }
                  }
                }
              }
            }
          }
        }
      }
    };
  public:
    AOInt(std::shared_ptr<const Geometry> g, const double gamma = 0.0, const bool ykw = false, const bool aux0 = false, const bool aux1 = false)
     : geom_(g), dim_(g->nbasis()), dim0_(aux0 ? g->naux() : g->nbasis()), dim1_(aux1 ? g->naux() : g->nbasis()), gamma_(gamma),
       aux0_(aux0), aux1_(aux1), yukawa_(ykw) {
      const size_t n = dim_;
      data_ = std::unique_ptr<double[]>(new double[n*n*dim0_*dim1_]);
      if (gamma != 0.0) {
        data2_ = std::unique_ptr<double[]>(new double[n*n*dim0_*dim1_]);
      }
      init();
    };
    ~AOInt() {};

    std::shared_ptr<File4> data()  { return std::shared_ptr<File4>(new File4(data_, dim_, dim0_, dim1_)); }; 
    std::shared_ptr<File4> data2() { return std::shared_ptr<File4>(new File4(data2_, dim_, dim0_, dim1_)); }; 

};

#endif


