//
// BAGEL - Parallel electron correlation program.
// Filename: tensor.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/tensor.h>
#include <src/smith/storagekramers.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template <typename DataType>
Tensor_<DataType>::Tensor_(vector<IndexRange> in, const bool kramers) : range_(in), rank_(in.size()), initialized_(false) {
  // make block list
  if (!in.empty()) {
    LoopGenerator lg(in);
    vector<vector<Index>> index = lg.block_loop();

    // first compute hashtags and length
    map<size_t, size_t> hashmap;
    size_t off = 0;
    for (auto& i : index) {
      size_t size = 1lu;
      vector<size_t> h;
      for (auto& j : i) {
        size *= j.size();
        h.push_back(j.key());
      }
      hashmap.emplace(generate_hash_key(h), size);
      off += size;
    }

    if (!kramers)
      data_ = make_shared<Storage<DataType>>(hashmap, false);
    else
      data_ = make_shared<StorageKramers<DataType>>(hashmap, false);
  } else {
    rank_ = 0;
    map<size_t, size_t> hashmap {{generate_hash_key(), 1lu}};
    data_ = make_shared<Storage<DataType>>(hashmap, false);
  }
}


template <typename DataType>
size_t Tensor_<DataType>::size_alloc() const {
  size_t out = 0lu;
  LoopGenerator lg(range_);
  vector<vector<Index>> index = lg.block_loop();
  for (auto& i : index) {
    vector<size_t> h;
    for (auto& j : i)
      h.push_back(j.key());
    out += data_->blocksize_alloc(h);
  }
  return out;
}


template <typename DataType>
vector<DataType> Tensor_<DataType>::diag() const {
  if (rank_ != 2 || range_.at(0) != range_.at(1))
    throw logic_error("Tensor_<DataType>::diag can be called only with a square tensor of rank 2");
  const size_t size = range_.at(0).back().offset() + range_.at(0).back().size();
  vector<DataType> buf(size);
  for (auto& i : range_.at(0)) {
    unique_ptr<DataType[]> data0 = get_block(i, i);
    for (int j = 0; j != i.size(); ++j) {
      buf[j+i.offset()] = data0[j+j*i.size()];
    }
  }
  return buf;
}


template <typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type> Tensor_<DataType>::matrix() const {
  vector<IndexRange> o = indexrange();
  assert(o.size() == 2);
  const int dim0 = o[0].size();
  const int dim1 = o[1].size();
  const int off0 = o[0].front().offset();
  const int off1 = o[1].front().offset();

  auto out = make_shared<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>(dim0, dim1);

  for (auto& i1 : o[1].range()) {
    for (auto& i0 : o[0].range()) {
      if (get_size_alloc(i0, i1)) {
        unique_ptr<DataType[]> target = get_block(i0, i1);
        out->copy_block(i0.offset()-off0, i1.offset()-off1, i0.size(), i1.size(), target.get());
      }
    }
  }
  return out;
}


// TODO parallelization
template <typename DataType>
shared_ptr<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type> Tensor_<DataType>::matrix2() const {
  const vector<IndexRange> o = indexrange();
  assert(o.size() == 4);

  const int dim0 = o[0].size();
  const int dim1 = o[1].size();
  const int dim2 = o[2].size();
  const int dim3 = o[3].size();
  const int off0 = o[0].front().offset();
  const int off1 = o[1].front().offset();
  const int off2 = o[2].front().offset();
  const int off3 = o[3].front().offset();

  auto out = make_shared<typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type>(dim0*dim1,dim2*dim3);
  for (auto& i3 : o[3].range()) {
    for (auto& i2 : o[2].range()) {
      for (auto& i1 : o[1].range()) {
        for (auto& i0 : o[0].range()) {
          if (get_size_alloc(i0, i1, i2, i3)) {
            unique_ptr<DataType[]> target = get_block(i0, i1, i2, i3);
            const DataType* ptr = target.get();
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ptr += i0.size())
                  copy_n(ptr, i0.size(), out->element_ptr(i0.offset()-off0+dim0*(j1-off1), j2-off2+dim2*(j3-off3)));
          }
        }
      }
    }
  }
  return out;
}


template <typename DataType>
shared_ptr<Civector<DataType>> Tensor_<DataType>::civec(shared_ptr<const Determinants> det) const {
  vector<IndexRange> o = indexrange();
  assert(o.size() == 1);

  int dim0 = 0;
  for (auto& i0 : o[0].range()) dim0 += i0.size();

  unique_ptr<DataType[]> civ(new DataType[dim0]);
  auto out = make_shared<Civector<DataType>>(det);
  for (auto& i0 : o[0].range()) {
    unique_ptr<DataType[]> target = get_block(i0);
    copy_n(target.get(), i0.size(), civ.get()+i0.offset());
  }

  copy_n(civ.get(), dim0, out->data());

  return out;
}


template <typename DataType>
void Tensor_<DataType>::print2(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 2);
  for (auto& i1 : o[1].range()) {
    for (auto& i0 : o[0].range()) {
      if (!this->get_size(i0, i1)) continue;
      unique_ptr<DataType[]> data = this->get_block(i0, i1);
      size_t iall = 0;
      for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
        for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
          if (abs(data[iall]) > thresh) {
            cout << "   " << setw(4) << j0 << " " << setw(4) << j1 <<
                           " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
          }
        }
    }
  }
  cout << "======================================" << endl << endl;
}


template <typename DataType>
void Tensor_<DataType>::print3(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 3);
  for (auto& i2 : o[2].range()) {
    for (auto& i1 : o[1].range()) {
      for (auto& i0 : o[0].range()) {
        if (!this->get_size(i0, i1, i2)) continue;
        unique_ptr<DataType[]> data = this->get_block(i0, i1, i2);
        size_t iall = 0;
        for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
              if (abs(data[iall]) > thresh) {
                cout << "   " << setw(4) << j0 << " " << setw(4) << j1 << " " << setw(4) << j2 <<
                               " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
              }
            }
      }
    }
  }
  cout << "======================================" << endl << endl;
}


template <typename DataType>
void Tensor_<DataType>::print4(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 4);
  for (auto& i3 : o[3].range()) {
    for (auto& i2 : o[2].range()) {
      for (auto& i1 : o[1].range()) {
        for (auto& i0 : o[0].range()) {
          if (!this->get_size_alloc(i0, i1, i2, i3)) continue;
          unique_ptr<DataType[]> data = this->get_block(i0, i1, i2, i3);
          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                  if (abs(data[iall]) > thresh) {
                     cout << "   " << setw(4) << j0 << " " << setw(4) << j1 <<
                                    " " << setw(4) << j2 << " " << setw(4) << j3 <<
                                    " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
                  }
                }
        }
      }
    }
  }
  cout << "======================================" << endl << endl;
}


template <typename DataType>
void Tensor_<DataType>::print5(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 5);
  for (auto& i4 : o[4].range()) {
    for (auto& i3 : o[3].range()) {
      for (auto& i2 : o[2].range()) {
        for (auto& i1 : o[1].range()) {
          for (auto& i0 : o[0].range()) {
            if (!this->get_size(i0, i1, i2, i3, i4)) continue;
            unique_ptr<DataType[]> data = this->get_block(i0, i1, i2, i3, i4);
            size_t iall = 0;
            for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
               for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                 for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                   for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                     for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                       if (abs(data[iall]) > thresh) {
                          cout << "   " << setw(4) << j0 << " " << setw(4) << j1 <<
                                         " " << setw(4) << j2 << " " << setw(4) << j3 << " " << setw(4) << j4 <<
                                         " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
                       }
                     }
          }
        }
      }
    }
  }
  cout << "======================================" << endl << endl;
}


template <typename DataType>
void Tensor_<DataType>::print6(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 6);
  for (auto& i5 : o[5].range()) {
    for (auto& i4 : o[4].range()) {
       for (auto& i3 : o[3].range()) {
         for (auto& i2 : o[2].range()) {
           for (auto& i1 : o[1].range()) {
             for (auto& i0 : o[0].range()) {
               if (!this->get_size(i0, i1, i2, i3, i4, i5)) continue;
               unique_ptr<DataType[]> data = this->get_block(i0, i1, i2, i3, i4, i5);
               size_t iall = 0;
               for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                 for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                   for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                     for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                       for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                         for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                           if (abs(data[iall]) > thresh) {
                              cout << "   " << setw(4) << j0 << " " << setw(4) << j1 <<
                                             " " << setw(4) << j2 << " " << setw(4) << j3 <<
                                             " " << setw(4) << j4 << " " << setw(4) << j5 <<
                                             " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
                           }
                         }
            }
          }
        }
      }
    }
  }
  cout << "======================================" << endl << endl;
}


template <typename DataType>
void Tensor_<DataType>::print8(string label, const double thresh) const {
  cout << endl << "======================================" << endl;
  cout << " > debug print out " << label << endl << endl;

  vector<IndexRange> o = indexrange();
  assert(o.size() == 8);
  for (auto& i7 : o[7].range()) {
    for (auto& i6 : o[6].range()) {
      for (auto& i5 : o[5].range()) {
        for (auto& i4 : o[4].range()) {
          for (auto& i3 : o[3].range()) {
            for (auto& i2 : o[2].range()) {
              for (auto& i1 : o[1].range()) {
                for (auto& i0 : o[0].range()) {
                  if (!this->get_size(i0, i1, i2, i3, i4, i5, i6, i7)) continue;
                  unique_ptr<DataType[]> data = this->get_block(i0, i1, i2, i3, i4, i5, i6, i7);
                  size_t iall = 0;
                  for (int j7 = i7.offset(); j7 != i7.offset()+i7.size(); ++j7)
                    for (int j6 = i6.offset(); j6 != i6.offset()+i6.size(); ++j6)
                      for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                        for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall) {
                                  if (abs(data[iall]) > thresh) {
                                     cout << "   " << setw(4) << j0 << " " << setw(4) << j1 <<
                                                    " " << setw(4) << j2 << " " << setw(4) << j3 <<
                                                    " " << setw(4) << j4 << " " << setw(4) << j5 <<
                                                    " " << setw(4) << j6 << " " << setw(4) << j7 <<
                                                    " " << setprecision(10) << setw(20) << fixed << data[iall] << endl;
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
  cout << "======================================" << endl << endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class Tensor_<double>;
template class Tensor_<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
