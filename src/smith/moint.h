//
// Newint - Parallel electron correlation program.
// Filename: moint.h
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

// A base class for electorn Correlation methods
// Certain kinds of MO integrals are formed.
//   - aaii (assumes DF - TODO half transformed DF vector might be available..)
//

#ifndef __SRC_SMITH_MOINT_H
#define __SRC_SMITH_MOINT_H

#include <memory>
#include <stdexcept>
#include <src/wfn/reference.h>
#include <src/smith/tensor.h>
#include <src/scf/fock.h>

namespace SMITH {

// the template parameter T specifies the storage type 

template <typename T>
class K2ext {
  protected:
    std::shared_ptr<Reference> ref_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor<T> > data_;

    // some handwritten drivers
    std::map<size_t, std::shared_ptr<DF_Full> > generate_list() {
      std::shared_ptr<DensityFit> df = ref_->geom()->df();
      std::shared_ptr<Coeff> coeff = ref_->coeff(); 

      // It is the easiest to do integral transformation for each blocks.
      assert(blocks_.size() == 4);
      std::map<size_t, std::shared_ptr<DF_Full> > dflist;
      // AO dimension
      const size_t nbasis = df->nbasis0();
      assert(df->nbasis0() == df->nbasis1());

      // TODO this part should be heavily parallelized
      // Also need to think a bit on the data layout. 
      // closed loop
      size_t cnt = blocks_[0].keyoffset();
      for (auto i0 = blocks_[0].range().begin(); i0 != blocks_[0].range().end(); ++i0, ++cnt) {
        std::shared_ptr<DF_Half> df_half = df->compute_half_transform(coeff->data()+nbasis*i0->offset(), i0->size()
                                                                      )->apply_J();
        // virtual loop
        size_t cnt2 = blocks_[1].keyoffset();
        for (auto i1 = blocks_[1].range().begin(); i1 != blocks_[1].range().end(); ++i1, ++cnt2) {
          std::shared_ptr<DF_Full> df_full = df_half->compute_second_transform(coeff->data()+nbasis*i1->offset(), i1->size());

          std::vector<size_t> h(1,cnt);
          h.push_back(cnt2);
          dflist.insert(make_pair(generate_hash_key(h), df_full));
        }
      }
      return dflist;
    };

    void form_4index(const std::map<size_t, std::shared_ptr<DF_Full> >& dflist) {
      // form four-index integrals
      // TODO this part should be heavily parallelized
      // TODO i01 < i23 symmetry should be used.
      size_t j0 = blocks_[0].keyoffset();
      for (auto i0 = blocks_[0].range().begin(); i0 != blocks_[0].range().end(); ++i0, ++j0) {
        size_t j1 = blocks_[1].keyoffset();
        for (auto i1 = blocks_[1].range().begin(); i1 != blocks_[1].range().end(); ++i1, ++j1) {
          // find three-index integrals
          std::vector<size_t> i01;
          i01.push_back(j0);
          i01.push_back(j1);
          auto iter01 = dflist.find(generate_hash_key(i01));
          assert(iter01 != dflist.end()); 
          std::shared_ptr<DF_Full> df01 = iter01->second; 
          size_t hashkey01 = generate_hash_key(i01);

          size_t j2 = blocks_[2].keyoffset();
          for (auto i2 = blocks_[2].range().begin(); i2 != blocks_[2].range().end(); ++i2, ++j2) {
            size_t j3 = blocks_[3].keyoffset();
            for (auto i3 = blocks_[3].range().begin(); i3 != blocks_[3].range().end(); ++i3, ++j3) {
              // find three-index integrals
              std::vector<size_t> i23;
              i23.push_back(j2);
              i23.push_back(j3);
              size_t hashkey23 = generate_hash_key(i23);
              if (hashkey23 > hashkey01) continue;

              auto iter23 = dflist.find(generate_hash_key(i23));
              assert(iter23 != dflist.end());
              std::shared_ptr<const DF_Full> df23 = iter23->second; 

              const size_t size = i0->size() * i1->size() * i2->size() * i3->size();

              // contract
              std::unique_ptr<double[]> target(new double[size]);
              df01->form_4index(target, df23); 

              // move in place
              std::vector<size_t> hash;
              hash.push_back(j0);
              hash.push_back(j1);
              hash.push_back(j2);
              hash.push_back(j3);

              if (hashkey23 != hashkey01) {
                std::unique_ptr<double[]> target2(new double[size]);
                const int s01 = i0->size() * i1->size();
                const int s23 = i2->size() * i3->size();
                mytranspose_(target.get(), &s01, &s23, target2.get()); 
                std::vector<size_t> hash2;
                hash2.push_back(j2);
                hash2.push_back(j3);
                hash2.push_back(j0);
                hash2.push_back(j1);
                data_->put_block(hash2, target2);
              }

              data_->put_block(hash, target);
            }
          }
        }
      }
    }; // vaaii_;

  public:
    K2ext(std::shared_ptr<Reference> r, std::vector<IndexRange> b) : ref_(r), blocks_(b) {
      // so far MOInt can be called for 2-external K integral and all-internals.
      if (blocks_[0] != blocks_[2] || blocks_[1] != blocks_[3]) 
        throw std::logic_error("MOInt called with wrong blocks");
      std::shared_ptr<Tensor<T> > tmp(new Tensor<T>(blocks_));
      data_ = tmp;
      form_4index(generate_list());
    };

    ~K2ext() {};

    std::shared_ptr<Tensor<T> > data() { return data_; };
    std::shared_ptr<Tensor<T> > tensor() { return data_; };

};


template <typename T>
class MOFock {
  protected:
    std::shared_ptr<Reference> ref_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor<T> > data_;
  public:
    MOFock(std::shared_ptr<Reference> r, std::vector<IndexRange> b) : ref_(r), blocks_(b) {
      // for simplicity, I assume that the Fock matrix is formed at once (may not be needed).
      assert(b.size() == 2 && b[0] == b[1]);

      std::shared_ptr<Tensor<T> > tmp(new Tensor<T>(blocks_));
      data_ = tmp;

      // TODO parallel not considered yet at all...
      std::shared_ptr<Fock<1> > fock0(new Fock<1>(ref_->geom(), ref_->hcore()));
      std::shared_ptr<Matrix1e> den = ref_->coeff()->form_density_rhf(ref_->nclosed()+ref_->nact());
      std::shared_ptr<Fock<1> > fock1(new Fock<1>(ref_->geom(), fock0, den, r->schwarz()));
      Matrix1e f = *r->coeff() % *fock1 * *r->coeff();
      size_t j0 = blocks_[0].keyoffset();
      for (auto i0 = blocks_[0].range().begin(); i0 != blocks_[0].range().end(); ++i0, ++j0) {
        size_t j1 = blocks_[1].keyoffset();
        for (auto i1 = blocks_[1].range().begin(); i1 != blocks_[1].range().end(); ++i1, ++j1) {
          const size_t size = i0->size() * i1->size();
          std::unique_ptr<double[]> target(new double[size]);
          double* buf = target.get();
          for (int k0 = i0->offset(); k0 != i0->offset()+i0->size(); ++k0) {
            for (int k1 = i1->offset(); k1 != i1->offset()+i1->size(); ++k1, ++buf) {
              *buf = f.element(k1, k0);
            }
          }

          std::vector<size_t> hash;
          hash.push_back(j1);
          hash.push_back(j0);
          data_->put_block(hash, target);
        }
      }
    };
    ~MOFock() {};

    std::shared_ptr<Tensor<T> > data() { return data_; };
    std::shared_ptr<Tensor<T> > tensor() { return data_; };
};

}

#endif

