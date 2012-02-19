//
// Author : Toru Shiozaki
// Date   : Feb 2012
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
      const size_t nbasis = df->nbasis();

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
          const size_t i1off = ref_->nclosed() + ref_->nact() + i1->offset();
          std::shared_ptr<DF_Full> df_full = df_half->compute_second_transform(coeff->data()+nbasis*i1off, i1->size());

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

          size_t j2 = blocks_[2].keyoffset();
          for (auto i2 = blocks_[2].range().begin(); i2 != blocks_[2].range().end(); ++i2, ++j2) {
            size_t j3 = blocks_[3].keyoffset();
            for (auto i3 = blocks_[3].range().begin(); i3 != blocks_[3].range().end(); ++i3, ++j3) {
              // find three-index integrals
              std::vector<size_t> i23;
              i23.push_back(j2);
              i23.push_back(j3);
              auto iter23 = dflist.find(generate_hash_key(i23));
              assert(iter23 != dflist.end());
              std::shared_ptr<DF_Full> df23 = iter23->second; 

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

};


template <typename T>
class MOFock {
  protected:
    std::shared_ptr<Reference> ref_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor<T> > data_;
  public:
    MOFock(std::shared_ptr<Reference> r, std::vector<IndexRange> b) : ref_(r), blocks_(b) {
    };
    ~MOFock() {};
};

}

#endif

