//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: df.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/df/df.h>
#include <src/df/dfdistt.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;
using namespace btas;


shared_ptr<DFDist> DFDist::copy() const {
  auto out = make_shared<DFDist>(df_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFDist> DFDist::clone() const {
  auto out = make_shared<DFDist>(df_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


tuple<int, vector<shared_ptr<const Shell>>> DFDist::get_ashell(const vector<shared_ptr<const Shell>>& all) {
  int out1;
  vector<shared_ptr<const Shell>> out2;
  // TODO without *2, H does not work. Perhaps need to think a bit more
  if (mpi__->size()*2 < all.size()) {
    int start, end;
    StaticDist d(naux_, mpi__->size());
    tie(start, end) = d.range(mpi__->rank());
    int num = 0;
    for (auto& s : all) {
      if (num >= start && num < end) {
        if (out2.empty()) out1 = num;
        out2.push_back(s);
      }
      num += s->nbasis();
    }
  } else {
    cout << endl << "   *** Warning *** Since the number of auxiliary shells is too small, we do not parallelize the Fock builder." << endl << endl;
    out1 = 0;
    out2 = all;
    serial_ = true;
  }

  return tie(out1, out2);
}


void DFDist::compute_2index(const vector<shared_ptr<const Shell>>& ashell, const double throverlap, const bool compute_inverse) {
  Timer time;

  // generates a task of integral evaluations
  TaskQueue<DFIntTask_OLD<DFDist>> tasks(ashell.size()*ashell.size());

  data2_ = make_shared<Matrix>(naux_, naux_, serial_);
  auto b3 = make_shared<const Shell>(ashell.front()->spherical());

  // naive static distribution
  int u = 0;
  int o0 = 0;
  for (auto& b0 : ashell) {
    int o1 = 0;
    for (auto& b1 : ashell) {
      if (o0 <= o1 && ((u++ % mpi__->size() == mpi__->rank()) || serial_))
        tasks.emplace_back(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, array<int,2>{{o0, o1}}, this);
      o1 += b1->nbasis();
    }
    o0 += b0->nbasis();
  }

  // these shell loops will be distributed across threads
  tasks.compute();

  if (!serial_)
    data2_->allreduce();

  time.tick_print("2-index ints");

  if (compute_inverse) {
    data2_->inverse_half(throverlap);
    // will use data2_ within node
    data2_->localize();
    time.tick_print("computing inverse");
  }
}


shared_ptr<const StaticDist> DFDist::make_table(const size_t astart) {
  vector<size_t> rec(mpi__->size());
  fill(rec.begin(), rec.end(), 0);

  mpi__->allgather(&astart, 1, rec.data(), 1);
  rec.push_back(naux_);

  return make_shared<const StaticDist>(rec);
}


pair<const double*, shared_ptr<RysInt>> DFDist::compute_batch(array<shared_ptr<const Shell>,4>& input) {
#ifdef LIBINT_INTERFACE
  shared_ptr<RysInt> eribatch = make_shared<Libint>(input);
#else
  shared_ptr<RysInt> eribatch = make_shared<ERIBatch>(input, 2.0);
#endif
  eribatch->compute();
  return {eribatch->data(), eribatch};
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform(const MatView c) const {
  const int nocc = c.extent(1);
  auto out = make_shared<DFHalfDist>(df_ ? df_ : shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_second(c));
  return out;
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform_swap(const MatView c) const {
  const int nocc = c.extent(1);
  auto out = make_shared<DFHalfDist>(df_ ? df_ : shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c)->swap());
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFHalfDist::compute_second_transform(const MatView c) const {
  const int nocc = c.extent(1);
  auto out = make_shared<DFFullDist>(df_, nindex1_, nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::copy() const {
  auto out = make_shared<DFHalfDist> (df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::clone() const {
  auto out = make_shared<DFHalfDist>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::merge_b1(shared_ptr<DFHalfDist> o) const {
  assert(naux() == o->naux() && nindex2() == o->nindex2());
  assert(block().size() == o->block().size());

  auto out = make_shared<DFHalfDist>(df_, nindex1() + o->nindex1());
  for (int i = 0; i != block_.size(); ++i)
    out->add_block(block(i)->merge_b1(o->block(i)));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::slice_b1(const int slice_start, const int slice_size) const {
  auto out = make_shared<DFHalfDist>(df_, slice_size);
  for (int i = 0; i != block_.size(); ++i)
    out->add_block(block(i)->slice_b1(slice_start, slice_size));
  return out;
}


shared_ptr<DFDist> DFHalfDist::back_transform(const MatView c) const{
  assert(df_->nindex1() == c.extent(0));
  auto out = make_shared<DFDist>(df_);
  for (auto& i : block_)
    out->add_block(i->transform_second(c, true));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::transform_occ(const shared_ptr<const Matrix> d) const {
  assert(nindex1_ == d->ndim());
  auto out = make_shared<DFHalfDist>(df_, d->mdim());
  for (auto& i : block_)
    out->add_block(i->transform_second(*d));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::apply_density(const shared_ptr<const Matrix> den) const {
  assert(den->mdim() == nindex2_);
  auto out = make_shared<DFHalfDist>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->transform_third(*den));
  return out;
}


shared_ptr<Matrix> DFHalfDist::compute_Kop_1occ(const shared_ptr<const Matrix> den, const double a) const {
  return apply_density(den)->form_2index(df_, a);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFFullDist::copy() const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::clone() const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::transform_occ1(const shared_ptr<const Matrix> d) const {
  assert(nindex1_ == d->ndim());
  auto out = make_shared<DFFullDist>(df_, d->mdim(), nindex2_);
  for (auto& i : block_)
    out->add_block(i->transform_second(*d));
  return out;
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DFHalfDist> DFFullDist::back_transform(const MatView c) const {
  assert(c.extent(0) == df_->nindex2());
  auto out = make_shared<DFHalfDist>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->transform_third(c, true));
  return out;
}


// 2RDM contractions
shared_ptr<DFFullDist> DFFullDist::apply_closed_2RDM(const double scale_exch) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_rhf_2RDM(scale_exch));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_uhf_2RDM(const btas::Tensor2<double>& amat, const btas::Tensor2<double>& bmat) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_uhf_2RDM(amat, bmat));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_2RDM(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const btas::Tensor4<double>& rdm) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_2RDM(rdm));
  return out;
}

shared_ptr<DFFullDist> DFFullDist::apply_2rdm_tran(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_2RDM_tr(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<Matrix> DFFullDist::form_aux_2index_apply_J(const shared_ptr<const DFFullDist> o, const double a) const {
  shared_ptr<Matrix> tmp = ParallelDF::form_aux_2index(o, a);
  return make_shared<Matrix>(*tmp * *df_->data2());
}


shared_ptr<Matrix> DFFullDist::form_4index_1fixed(const shared_ptr<const DFFullDist> o, const double a, const size_t n) const {
  // TODO needs more work
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = block_[0]->form_4index_1fixed(o->block_[0], a, n);
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<Matrix> DFFullDist::form_4index_diagonal() const {
  // TODO needs more work
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = block_[0]->form_4index_diagonal();
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<Matrix> DFFullDist::form_4index_diagonal_part() const {
  // TODO needs more work
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = block_[0]->form_4index_diagonal_part();
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<DFFullDist> DFFullDist::swap() const {
  auto out = make_shared<DFFullDist>(shared_from_this(), nocc2(), nocc1());
  for (auto& i : block_)
    out->add_block(i->swap());
  return out;
}


//////// apply J functions ////////

shared_ptr<DFFullDist> DFFullDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFFullDist> out = clone();
#ifdef HAVE_MPI_H
  if (!serial_) {
    Timer mult(3);
    auto work = make_shared<DFDistT>(shared_from_this());
    mult.tick_print("Form DFDistT");
    work = work->apply_J(d);
    mult.tick_print("Application of Inverse");
    work->get_paralleldf(out);
    mult.tick_print("Return DFDist");
  } else {
#else
  {
#endif
    auto j = block_.begin();
    for (auto& i : out->block_) {
      i->zero();
      i->contrib_apply_J(*j, d);
      ++j;
    }
  }
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFHalfDist> out = clone();
#ifdef HAVE_MPI_H
  if (!serial_) {
    Timer mult(3);
    auto work = make_shared<DFDistT>(shared_from_this());
    mult.tick_print("Form DFDistT");
    work = work->apply_J(d);
    mult.tick_print("Application of Inverse");
    work->get_paralleldf(out);
    mult.tick_print("Return DFDist");
  } else {
#else
  {
#endif
    auto j = block_.begin();
    for (auto& i : out->block_) {
      i->zero();
      i->contrib_apply_J(*j, d);
      ++j;
    }
  }
  return out;
}
