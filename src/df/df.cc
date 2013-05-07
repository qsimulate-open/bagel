//
// BAGEL - Parallel electron correlation program.
// Filename: df.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#include <src/df/df.h>
#include <src/df/dfdistt.h>
#include <src/rysint/eribatch.h>
#include <src/rysint/libint.h>

using namespace std;
using namespace bagel;


ParallelDF::ParallelDF(const size_t naux, const size_t nb1, const size_t nb2, std::shared_ptr<const ParallelDF> df, std::shared_ptr<Matrix> dat)
 : naux_(naux), nindex1_(nb1), nindex2_(nb2), df_(df), data2_(dat), serial_(df ? df->serial_ : false) {

}


shared_ptr<Matrix> ParallelDF::form_2index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = (!swap) ? block_[0]->form_2index(o->block_[0], a) : o->block_[0]->form_2index(block_[0], a);
  if (!serial_)
    out->allreduce();
  return out;
}


unique_ptr<double[]> ParallelDF::form_4index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  unique_ptr<double[]> out = (!swap) ? block_[0]->form_4index(o->block_[0], a) : o->block_[0]->form_4index(block_[0], a);

  // all reduce
  const size_t size = block_[0]->b2size()*o->block_[0]->b2size() * block_[0]->b1size()*o->block_[0]->b1size();
  if (!serial_)
    mpi__->allreduce(out.get(), size);
  return out;
}


shared_ptr<Matrix> ParallelDF::form_aux_2index(shared_ptr<const ParallelDF> o, const double a) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
#ifdef HAVE_MPI_H
  if (!serial_) {
    auto work = make_shared<DFDistT>(shared_from_this());
    auto work2 = make_shared<DFDistT>(o);
    return work->form_aux_2index(work2, a).front();
  } else {
#else
  {
#endif
    return block_[0]->form_aux_2index(o->block_[0], a);
  }
}


void ParallelDF::daxpy(const double a, const shared_ptr<const ParallelDF> o) {
  assert(block_.size() == o->block_.size());
  auto j = o->block_.begin();
  for (auto& i : block_)
    i->daxpy(a, *j++);
}


void ParallelDF::scale(const double a) {
  for (auto& i : block_)
    i->scale(a);
}


void ParallelDF::add_block(shared_ptr<DFBlock> o) {
  block_.push_back(o);
}


unique_ptr<double[]> ParallelDF::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  // first thing is to find the node
  tuple<size_t, size_t> info = adist_now()->locate(i);

  // ask for the data to inode
  if (get<0>(info) == mpi__->rank()) {
    return block_[0]->get_block(i, id, j, jd, k, kd);
  } else {
    throw logic_error("ParallelDF::get_block is an intra-node function (or bug?)");
  }
  return unique_ptr<double[]>();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DFDist::add_direct_product(const vector<const double*> cd, const vector<const double*> dd, const double a) {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  if (cd.size() != dd.size()) throw logic_error("Illegal call of DFDist::DFDist");

  for (auto c = cd.begin(), d = dd.begin(); c != cd.end(); ++c, ++d)
    block_[0]->add_direct_product(*c+block_[0]->astart(), *d, a);
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
    for (auto iter = all.begin(); iter != all.end(); ++iter) {
      if (num >= start && num < end) {
        if (out2.empty()) out1 = num;
        out2.push_back(*iter);
      }
      num += (*iter)->nbasis();
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
  vector<DFIntTask_OLD<DFDist>> tasks;
  tasks.reserve(ashell.size()*ashell.size());

  data2_ = make_shared<Matrix>(naux_, naux_, serial_);
  auto b3 = make_shared<const Shell>(ashell.front()->spherical());

  // naive static distribution
  int u = 0;
  int o0 = 0;
  for (auto& b0 : ashell) {
    int o1 = 0;
    for (auto& b1 : ashell) {
      if (o0 <= o1 && ((u++ % mpi__->size() == mpi__->rank()) || serial_))
        tasks.push_back(DFIntTask_OLD<DFDist>(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, array<int,2>{{o0, o1}}, this));
      o1 += b1->nbasis();
    }
    o0 += b0->nbasis();
  }

  // these shell loops will be distributed across threads
  TaskQueue<DFIntTask_OLD<DFDist>> tq(tasks);
  tq.compute(resources__->max_num_threads());

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

  mpi__->allgather(&astart, 1, &rec[0], 1); 
  rec.push_back(naux_);

  return make_shared<const StaticDist>(rec);
}


pair<const double*, shared_ptr<RysInt>> DFDist::compute_batch(array<shared_ptr<const Shell>,4>& input) {
#ifdef LIBINT_INTERFACE
  auto eribatch = make_shared<Libint>(input);
#else
  auto eribatch = make_shared<ERIBatch>(input, 2.0);
#endif
  eribatch->compute();
  return make_pair(eribatch->data(), eribatch);
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const Matrix> den) const {
  return compute_Jop(shared_from_this(), den);
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const ParallelDF> o, const shared_ptr<const Matrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const Matrix> tmp0 = o->compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return compute_Jop_from_cd(tmp0);
}


shared_ptr<Matrix> ParallelDF::compute_Jop_from_cd(shared_ptr<const Matrix> tmp0) const {
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = block_[0]->form_mat(tmp0->data()+block_[0]->astart());
  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<Matrix> ParallelDF::compute_cd(const shared_ptr<const Matrix> den, shared_ptr<const Matrix> dat2, const bool onlyonce) const {
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = make_shared<Matrix>(naux_, 1, true);

  // D = (D|rs)*d_rs
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  unique_ptr<double[]> tmp = block_[0]->form_vec(den);
  copy_n(tmp.get(), block_[0]->asize(), tmp0->data()+block_[0]->astart());
  // All reduce
  if (!serial_)
    tmp0->allreduce();

  tmp0 = make_shared<Matrix>(*dat2 * *tmp0);
  if (!onlyonce)
    tmp0 = make_shared<Matrix>(*dat2 * *tmp0);
  return tmp0;
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform(const double* c, const size_t nocc) const {
  auto out = make_shared<DFHalfDist>(shared_from_this(), nocc);
  for (auto& i : block_) 
    out->add_block(i->transform_second(c, nocc));
  return out;
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform_swap(const double* c, const size_t nocc) const {
  auto out = make_shared<DFHalfDist>(shared_from_this(), nocc);
  for (auto& i : block_) 
    out->add_block(i->transform_third(c, nocc)->swap());
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFHalfDist::compute_second_transform(const double* c, const size_t nocc) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nocc);
  for (auto& i : block_) 
    out->add_block(i->transform_third(c, nocc));
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


shared_ptr<DFDist> DFHalfDist::back_transform(const double* c) const{
  auto out = make_shared<DFDist>(df_);
  for (auto& i : block_) 
    out->add_block(i->transform_second(c, df_->nindex1(), true));
  return out;
}


void DFHalfDist::rotate_occ(const double* d) {
  for (auto& i : block_)
    i = i->transform_second(d, nindex1_);
}


shared_ptr<DFHalfDist> DFHalfDist::apply_density(const double* den) const {
  auto out = make_shared<DFHalfDist>(df_, nindex1_); 
  for (auto& i : block_) 
    out->add_block(i->transform_third(den, nindex2_));
  return out;
}


shared_ptr<Matrix> DFHalfDist::compute_Kop_1occ(const double* den, const double a) const {
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


void DFFullDist::symmetrize() {
  for (auto& i : block_)
    i->symmetrize();
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DFHalfDist> DFFullDist::back_transform(const double* c) const {
  auto out = make_shared<DFHalfDist>(df_, nindex1_);
  for (auto& i : block_) 
    out->add_block(i->transform_third(c, df_->nindex2(), true));
  return out;
}


// 2RDM contractions
shared_ptr<DFFullDist> DFFullDist::apply_closed_2RDM(const double scale_exch) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_) 
    out->add_block(i->apply_rhf_2RDM(scale_exch));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_) 
    out->add_block(i->apply_uhf_2RDM(amat, bmat));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_) 
    out->add_block(i->apply_2RDM(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm) const {
  auto out = make_shared<DFFullDist>(df_, nindex1_, nindex2_);
  for (auto& i : block_) 
    out->add_block(i->apply_2RDM(rdm));
  return out;
}


shared_ptr<Matrix> DFFullDist::form_aux_2index_apply_J(const shared_ptr<const DFFullDist> o, const double a) const {
  shared_ptr<Matrix> tmp = ParallelDF::form_aux_2index(o, a);
  return make_shared<Matrix>(*tmp * *df_->data2());
}


unique_ptr<double[]> DFFullDist::form_4index_1fixed(const shared_ptr<const DFFullDist> o, const double a, const size_t n) const {
  // TODO needs more work
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  unique_ptr<double[]> out = block_[0]->form_4index_1fixed(o->block_[0], a, n);
  const size_t size = block_[0]->b1size() * block_[0]->b2size() * o->block_[0]->b1size();
  if (!serial_)
    mpi__->allreduce(out.get(), size);
  return out;
}


void DFFullDist::set_product(const shared_ptr<const DFFullDist> o, const unique_ptr<double[]>& c, const int jdim, const size_t off) {
  // TODO needs more work
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  block_[0]->copy_block(o->block_[0]->form_Dj(c, jdim), jdim, off*block_[0]->asize());
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
    mult.tick_print("Applicatoin of Inverse");
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


