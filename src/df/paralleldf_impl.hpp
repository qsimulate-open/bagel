//
// BAGEL - Parallel electron correlation program.
// Filename: paralleldf_impl.hpp
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#ifdef PARALLELDF_HEADERS

#ifndef __SRC_DF_PARALLELDF_IMPL_HPP
#define __SRC_DF_PARALLELDF_IMPL_HPP


namespace bagel {

template <typename DataType, typename MatrixType, typename BlockType>
ParallelDFit<DataType,MatrixType,BlockType>::ParallelDFit(const size_t naux, const size_t nb1, const size_t nb2, std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> df, std::shared_ptr<MatrixType> dat)
 : naux_(naux), nindex1_(nb1), nindex2_(nb2), df_(df), data2_(dat), serial_(df ? df->serial_ : false) {

}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::form_2index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw std::logic_error("so far assumes block_.size() == 1");
  std::shared_ptr<MatrixType> out = (!swap) ? block_[0]->form_2index(o->block_[0], a) : o->block_[0]->form_2index(block_[0], a);
  if (!serial_)
    out->allreduce();
  return out;
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::form_4index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw std::logic_error("so far assumes block_.size() == 1");
  std::shared_ptr<MatrixType> out = (!swap) ? block_[0]->form_4index(o->block_[0], a) : o->block_[0]->form_4index(block_[0], a);

  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::form_aux_2index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw std::logic_error("so far assumes block_.size() == 1");
#ifdef HAVE_MPI_H
  if (!serial_) {
    auto work = std::make_shared<DFDistT>(this->shared_from_this());
    auto work2 = std::make_shared<DFDistT>(o);
    return work->form_aux_2index(work2, a).front();
  } else {
#else
  {
#endif
    return block_[0]->form_aux_2index(o->block_[0], a);
  }
}


template <typename DataType, typename MatrixType, typename BlockType>
void ParallelDFit<DataType,MatrixType,BlockType>::ax_plus_y(const double a, const std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o) {
  assert(block_.size() == o->block_.size());
  auto j = o->block_.begin();
  for (auto& i : block_)
    i->ax_plus_y(a, *j++);
}


template <typename DataType, typename MatrixType, typename BlockType>
void ParallelDFit<DataType,MatrixType,BlockType>::scale(const double a) {
  for (auto& i : block_)
    i->scale(a);
}


template <typename DataType, typename MatrixType, typename BlockType>
void ParallelDFit<DataType,MatrixType,BlockType>::symmetrize() {
  for (auto& i : block_)
    i->symmetrize();
}


template <typename DataType, typename MatrixType, typename BlockType>
void ParallelDFit<DataType,MatrixType,BlockType>::add_block(std::shared_ptr<BlockType> o) {
  block_.push_back(o);
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  if (block_.size() != 1) throw std::logic_error("so far assumes block_.size() == 1");
  // first thing is to find the node
  std::tuple<size_t, size_t> info = adist_now()->locate(i);

  // date has to be localised in this node
  if (std::get<0>(info) == mpi__->rank() && !block_[0]->averaged()) {
    return block_[0]->get_block(i, id, j, jd, k, kd);
  } else {
    throw std::logic_error("ParallelDFit<DataType,MatrixType,BlockType>::get_block is an intra-node function (or bug?)");
  }
  return nullptr;
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::compute_Jop_from_cd(std::shared_ptr<const MatrixType> tmp0) const {
  if (block_.size() != 1) throw std::logic_error("compute_Jop so far assumes block_.size() == 1");
  std::shared_ptr<MatrixType> out = block_[0]->form_mat(tmp0->data()+block_[0]->astart());
  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::compute_cd(const std::shared_ptr<const MatrixType> den, std::shared_ptr<const MatrixType> dat2, const bool onlyonce) const {
  if (!dat2 && !data2_) throw std::logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = std::make_shared<MatrixType>(naux_, 1, true);

  // D = (D|rs)*d_rs
  if (block_.size() != 1) throw std::logic_error("compute_Jop so far assumes block_.size() == 1");
  std::unique_ptr<DataType[]> tmp = block_[0]->form_vec(den);
  std::copy_n(tmp.get(), block_[0]->asize(), tmp0->data()+block_[0]->astart());
  // All reduce
  if (!serial_)
    tmp0->allreduce();

  tmp0 = std::make_shared<MatrixType>(*dat2 * *tmp0);
  if (!onlyonce)
    tmp0 = std::make_shared<MatrixType>(*dat2 * *tmp0);
  return tmp0;
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::compute_Jop(const std::shared_ptr<const MatrixType> den) const {
  return compute_Jop(this->shared_from_this(), den);
}


template <typename DataType, typename MatrixType, typename BlockType>
std::shared_ptr<MatrixType> ParallelDFit<DataType,MatrixType,BlockType>::compute_Jop(const std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const std::shared_ptr<const MatrixType> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  std::shared_ptr<const MatrixType> tmp0 = o->compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return compute_Jop_from_cd(tmp0);
}


}


#endif
#endif
