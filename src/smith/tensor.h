//
// BAGEL - Parallel electron correlation program.
// Filename: tensor.h
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


#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <stddef.h>
#include <list>
#include <map>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <src/fci/civec.h>
#include <type_traits>
#include <src/math/matrix.h>
#include <src/math/matop.h>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>

namespace bagel {
namespace SMITH {

// this assumes < 256 blocks; TODO runtime determination?
const static int shift = 8;

/* obsolete function */
static
size_t generate_hash_key(const std::vector<size_t>& o) {
  size_t out = 0;
  for (auto i = o.rbegin(); i != o.rend(); ++i) {
    out <<= shift;
    out += *i;
  }
  return out;
}

static
size_t generate_hash_key() { return 0; }

template<class T, typename... args>
size_t generate_hash_key(const T& head, const args&... tail) {
  return (generate_hash_key(tail...) << shift) + head.key();
}


template <typename T>
class Tensor {
  protected:

    std::vector<IndexRange> range_;
    std::shared_ptr<T> data_;
    int rank_;
    bool initialized_;

  public:
    Tensor(std::vector<IndexRange> in, bool init = true) : range_(in), rank_(in.size()), initialized_(init) {
      // make blocl list
      if (!in.empty()) {
        LoopGenerator lg(in);
        std::vector<std::vector<Index>> index = lg.block_loop();

        // first compute hashtags and length
        std::map<size_t, size_t> hashmap;
        size_t off = 0;
        for (auto& i : index) {
          size_t size = 1lu;
          std::vector<size_t> h;
          for (auto& j : i) {
            size *= j.size();
            h.push_back(j.key());
          }
          hashmap.insert(std::make_pair(generate_hash_key(h), size));
          off += size;
        }

        data_ = std::make_shared<T>(hashmap, init);
      } else {
        rank_ = 0;
        std::map<size_t, size_t> hashmap {{0lu, 1lu}};
        data_ = std::make_shared<T>(hashmap, init);
      }
    }

    void initialize() { data_->initialize(); }

    Tensor<T>& operator=(const Tensor<T>& o) {
      *data_ = *(o.data_);
      return *this;
    }

    std::shared_ptr<Tensor<T>> clone() const {
      return std::make_shared<Tensor<T>>(range_);
    }

    std::shared_ptr<Tensor<T>> copy() const {
      std::shared_ptr<Tensor<T>> out = clone();
      *out = *this;
      return out;
    }

    void ax_plus_y(const double a, const Tensor<T>& o) { data_->ax_plus_y(a, o.data_); }
    void ax_plus_y(const double a, const std::shared_ptr<Tensor<T>> o) { data_->ax_plus_y(a, o->data_); }

    void scale(const double a) { data_->scale(a); }

    double dot_product(const Tensor<T>& o) { return data_->dot_product(*o.data_); }
    double dot_product(const std::shared_ptr<Tensor<T>>& o) { return data_->dot_product(*o->data_); }

    size_t size() const { return data_->length(); }
    size_t length() const { return data_->length(); }

    double norm() { return std::sqrt(dot_product(*this)); }
    double rms() { return std::sqrt(dot_product(*this)/size()); }

    std::vector<IndexRange> indexrange() const { return range_; }

    template<typename ...args>
    std::unique_ptr<double[]> get_block(const args& ...p) const {
      return data_->get_block(generate_hash_key(p...));
    }

    template<typename ...args>
    std::unique_ptr<double[]> move_block(const args& ...p) const {
      return data_->move_block(generate_hash_key(p...));
    }

    template<typename ...args>
    void put_block(std::unique_ptr<double[]>& o, const args& ...p) const {
      data_->put_block(generate_hash_key(p...), o);
    }

    template<typename ...args>
    void add_block(std::unique_ptr<double[]>& o, const args& ...p) const {
      data_->add_block(generate_hash_key(p...), o);
    }

    template<typename ...args>
    size_t get_size(const args& ...p) const {
      return data_->blocksize(generate_hash_key(p...));
    }

/****************** following functions are obsolete *************************/
    std::unique_ptr<double[]> get_block(const std::vector<size_t>& p) const {
      assert(p.size() == rank_ || (rank_ == 0 && p.size() == 1));
      if (data_ == nullptr) throw std::logic_error("Tensor not initialized");
      return data_->get_block(generate_hash_key(p));
    }

    std::unique_ptr<double[]> get_block(const std::initializer_list<size_t>& p) const {
      return get_block(std::vector<size_t>(p.begin(), p.end()));
    }

    std::unique_ptr<double[]> move_block(const std::vector<size_t>& p) {
      assert(p.size() == rank_ || (rank_ == 0 && p.size() == 1));
      return data_->move_block(generate_hash_key(p));
    }

    std::unique_ptr<double[]> move_block(const std::initializer_list<size_t>& p) const {
      return move_block(std::vector<size_t>(p.begin(), p.end()));
    }

    void put_block(const std::vector<size_t>& p, std::unique_ptr<double[]>& o) {
      data_->put_block(generate_hash_key(p), o);
    }

    void put_block(const std::initializer_list<size_t>& p, std::unique_ptr<double[]>& o) {
      put_block(std::vector<size_t>(p.begin(), p.end()), o);
    }

    void add_block(const std::vector<size_t>& p, const std::unique_ptr<double[]>& o) {
      if (data_ == nullptr) throw std::logic_error("Tensor not initialized");
      data_->add_block(generate_hash_key(p), o);
    }

    void add_block(const std::initializer_list<size_t>& p, std::unique_ptr<double[]>& o) {
      add_block(std::vector<size_t>(p.begin(), p.end()), o);
    }

    size_t get_size(const std::vector<size_t>& p) {
      assert(p.size() == rank_);
      return data_->blocksize(generate_hash_key(p));
    }
/****************** to here *************************/

    void zero() {
      data_->zero();
    }

    std::vector<double> diag() {
      if (rank_ != 2 || range_[0] != range_[1])
        throw std::logic_error("Tensor::diag can be called only with a square tensor of rank 2");
      const size_t size = range_[0].back().offset() + range_[0].back().size();
      std::vector<double> buf(size);
      for (auto& i : range_[0]) {
        std::unique_ptr<double[]> data0 = move_block(i, i);
        for (int j = 0; j != i.size(); ++j) {
          buf[j+i.offset()] = data0[j+j*i.size()];
        }
        put_block(data0, i, i);
      }
      return buf;
    }


    std::shared_ptr<Tensor<T>> add_dagger() {
      std::shared_ptr<Tensor<T>> out = clone();
      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 4);
      for (auto& i3 : o[3].range()) {
        for (auto& i2 : o[2].range()) {
          for (auto& i1 : o[1].range()) {
            for (auto& i0 : o[0].range()) {
              std::unique_ptr<double[]>       data0 = get_block(i0, i1, i2, i3);
              const std::unique_ptr<double[]> data1 = get_block(i2, i3, i0, i1);
              sort_indices<2,3,0,1,1,1,1,1>(data1, data0, i2.size(), i3.size(), i0.size(), i1.size());
              out->put_block(data0, i0, i1, i2, i3);
            }
          }
        }
      }
      return out;
    }


    std::shared_ptr<Matrix> matrix() const {
      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 2);
      const int dim0 = o[0].size();
      const int dim1 = o[1].size();
      const int off0 = o[0].front().offset();
      const int off1 = o[1].front().offset();

      auto out = std::make_shared<Matrix>(dim0, dim1);

      for (auto& i1 : o[1].range()) {
        for (auto& i0 : o[0].range()) {
          std::unique_ptr<double[]> target = get_block(i0, i1);
          out->copy_block(i0.offset()-off0, i1.offset()-off1, i0.size(), i1.size(), target);
        }
      }
      return out;
    }


    // TODO parallelization
    std::shared_ptr<Matrix> matrix2() const {
      const std::vector<IndexRange> o = indexrange();
      assert(o.size() == 4);

      const int dim0 = o[0].size();
      const int dim1 = o[1].size();
      const int dim2 = o[2].size();
      const int dim3 = o[3].size();
      const int off0 = o[0].front().offset();
      const int off1 = o[1].front().offset();
      const int off2 = o[2].front().offset();
      const int off3 = o[3].front().offset();

      auto out = std::make_shared<Matrix>(dim0*dim1,dim2*dim3);
      for (auto& i3 : o[3].range()) {
        for (auto& i2 : o[2].range()) {
          for (auto& i1 : o[1].range()) {
            for (auto& i0 : o[0].range()) {
              std::unique_ptr<double[]> target = get_block(i0, i1, i2, i3);
              const double* ptr = target.get();
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ptr += i0.size())
                    std::copy_n(ptr, i0.size(), out->element_ptr(i0.offset()-off0+dim0*(j1-off1), j2-off2+dim2*(j3-off3)));
            }
          }
        }
      }
      return out;
    }


    std::shared_ptr<Civec> civec(std::shared_ptr<const Determinants> det) const {
      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 1);

      int dim0 = 0;
      for (auto& i0 : o[0].range()) dim0 += i0.size();

      std::unique_ptr<double[]> civ(new double[dim0]);
      auto out = std::make_shared<Civec>(det);
      for (auto& i0 : o[0].range()) {
        std::unique_ptr<double[]> target = get_block(i0);
        std::copy_n(target.get(), i0.size(), civ.get()+i0.offset());
      }

      std::copy_n(civ.get(), dim0, out->data());

      return out;
    }


    void print1(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 1);
      for (auto& i0 : o[0].range()) {
        if (!this->get_size(i0)) continue;
        std::unique_ptr<double[]> data = this->get_block(i0);
        size_t iall = 0;
        for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
          if (fabs(data[iall]) > thresh) {
            std::cout << "   " << std::setw(4) << j0 << " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }


    void print2(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 2);
      for (auto& i1 : o[1].range()) {
        for (auto& i0 : o[0].range()) {
          // if this block is not included in the current wave function, skip it
          if (!this->get_size(i0, i1)) continue;
          std::unique_ptr<double[]> data = this->get_block(i0, i1);
          size_t iall = 0;
          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
              if (fabs(data[iall]) > thresh) {
                std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 <<
                               " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
              }
            }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }

    void print3(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 3);
      for (auto& i2 : o[2].range()) {
        for (auto& i1 : o[1].range()) {
          for (auto& i0 : o[0].range()) {
            // if this block is not included in the current wave function, skip it
            if (!this->get_size(i0, i1, i2)) continue;
            std::unique_ptr<double[]> data = this->get_block(i0, i1, i2);
            size_t iall = 0;
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                  if (fabs(data[iall]) > thresh) {
                    std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 << " " << std::setw(4) << j2 <<
                                   " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
                  }
                }
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }


    void print4(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 4);
      for (auto& i3 : o[3].range()) {
        for (auto& i2 : o[2].range()) {
          for (auto& i1 : o[1].range()) {
            for (auto& i0 : o[0].range()) {
              // if this block is not included in the current wave function, skip it
              if (!this->get_size(i0, i1, i2, i3)) continue;
              std::unique_ptr<double[]> data = this->get_block(i0, i1, i2, i3);
              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                      if (fabs(data[iall]) > thresh) {
                         std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 <<
                                        " " << std::setw(4) << j2 << " " << std::setw(4) << j3 <<
                                        " " << std::setprecision(10) << std::setw(15) << std::fixed << data[(j0-i0.offset())+i0.size()*((j1-i1.offset())+i1.size()*((j2-i2.offset())+i2.size()*((j3-i3.offset()))))] << std::endl;  // testing
                                        //" " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
                      }
                    }
            }
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }

    void print5(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 5);
      for (auto& i4 : o[4].range()) {
        for (auto& i3 : o[3].range()) {
          for (auto& i2 : o[2].range()) {
            for (auto& i1 : o[1].range()) {
              for (auto& i0 : o[0].range()) {
                // if this block is not included in the current wave function, skip it
                if (!this->get_size(i0, i1, i2, i3, i4)) continue;
                std::unique_ptr<double[]> data = this->get_block(i0, i1, i2, i3, i4);
                size_t iall = 0;
                for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                   for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                     for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                       for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                         for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                           if (fabs(data[iall]) > thresh) {
                              std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 <<
                                             " " << std::setw(4) << j2 << " " << std::setw(4) << j3 << " " << std::setw(4) << j4 <<
                                             " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
                           }
                         }
              }
            }
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }

    void print6(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 6);
      for (auto& i5 : o[5].range()) {
        for (auto& i4 : o[4].range()) {
           for (auto& i3 : o[3].range()) {
             for (auto& i2 : o[2].range()) {
               for (auto& i1 : o[1].range()) {
                 for (auto& i0 : o[0].range()) {
                   // if this block is not included in the current wave function, skip it
                   if (!this->get_size(i0, i1, i2, i3, i4, i5)) continue;
                   std::unique_ptr<double[]> data = this->get_block(i0, i1, i2, i3, i4, i5);
                   size_t iall = 0;
                   for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                     for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                       for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                         for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                           for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                             for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                               if (fabs(data[iall]) > thresh) {
                                  std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 <<
                                                 " " << std::setw(4) << j2 << " " << std::setw(4) << j3 <<
                                                 " " << std::setw(4) << j4 << " " << std::setw(4) << j5 <<
                                                 " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
                               }
                             }
                }
              }
            }
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }

    void print8(std::string label, const double thresh = 5.0e-2) {
      std::cout << std::endl << "======================================" << std::endl;
      std::cout << " > debug print out " << label << std::endl << std::endl;

      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 8);
      for (auto& i7 : o[7].range()) {
        for (auto& i6 : o[6].range()) {
          for (auto& i5 : o[5].range()) {
            for (auto& i4 : o[4].range()) {
              for (auto& i3 : o[3].range()) {
                for (auto& i2 : o[2].range()) {
                  for (auto& i1 : o[1].range()) {
                    for (auto& i0 : o[0].range()) {
                      // if this block is not included in the current wave function, skip it
                      if (!this->get_size(i0, i1, i2, i3, i4, i5, i6, i7)) continue;
                      std::unique_ptr<double[]> data = this->get_block(i0, i1, i2, i3, i4, i5, i6, i7);
                      size_t iall = 0;
                      for (int j7 = i7.offset(); j7 != i7.offset()+i7.size(); ++j7)
                        for (int j6 = i6.offset(); j6 != i6.offset()+i6.size(); ++j6)
                          for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                            for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                                      if (fabs(data[iall]) > thresh) {
                                         std::cout << "   " << std::setw(4) << j0 << " " << std::setw(4) << j1 <<
                                                        " " << std::setw(4) << j2 << " " << std::setw(4) << j3 <<
                                                        " " << std::setw(4) << j4 << " " << std::setw(4) << j5 <<
                                                        " " << std::setw(4) << j6 << " " << std::setw(4) << j7 <<
                                                        " " << std::setprecision(10) << std::setw(15) << std::fixed << data[iall] << std::endl;
                                      }
                                    }
                    }
                  }
                }
              }
            }
          }
        }
      }
      std::cout << "======================================" << std::endl << std::endl;
    }
};

}
}

#endif

