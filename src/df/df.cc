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


#include <memory>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <list>
#include <src/rysint/eribatch.h>
#include <src/rysint/libint.h>
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/util/timer.h>
#include <src/df/df.h>
#include <src/parallel/paramatrix.h>
#include <src/parallel/mpi_interface.h>

#include <src/df/dfinttask_old.h>
#include <src/df/dfdistt.h>

using namespace std;
using namespace bagel;


ParallelDF::ParallelDF(const size_t naux, const size_t nb1, const size_t nb2)
 : naux_(naux), nindex1_(nb1), nindex2_(nb2) {

}


int ParallelDF::get_node(const int off) const {
  int out = df_->global_table_.upper_bound(off)->second;
  return out;
}


shared_ptr<Matrix> ParallelDF::form_2index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  shared_ptr<ParaMatrix> out = (!swap) ? block_->form_2index(o->block_, a) : o->block_->form_2index(block_, a);
  out->allreduce();
  return out;
}


unique_ptr<double[]> ParallelDF::form_4index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  unique_ptr<double[]> out = (!swap) ? block_->form_4index(o->block_, a) : o->block_->form_4index(block_, a);

  // all reduce
  const size_t size = block_->b2size()*o->block_->b2size() * block_->b1size()*o->block_->b1size();
  mpi__->allreduce(out.get(), size);
  return out;
}


shared_ptr<Matrix> ParallelDF::form_aux_2index(shared_ptr<const ParallelDF> o, const double a) const {
#ifdef HAVE_MPI_H
  shared_ptr<ParaMatrix> out(new ParaMatrix(naux_, naux_));
  shared_ptr<DFDistT> work(new DFDistT(shared_from_this()));
  shared_ptr<DFDistT> work2 = this == o.get() ? work : shared_ptr<DFDistT>(new DFDistT(o));
  assert(work->size() == work2->size());
  dgemm_("T", "N", naux_, naux_, work->size(), a, work->data(), work->size(), work2->data(), work2->size(), 0.0, out->data(), naux_); 
  out->allreduce();
  return out;
#else
  return block_->form_aux_2index(o->block_, a);
#endif
}


void ParallelDF::daxpy(const double a, const shared_ptr<const ParallelDF> o) {
  block_->daxpy(a, o->block_);
}


void ParallelDF::scale(const double a) {
  block_->scale(a);
}


void ParallelDF::add_block(shared_ptr<DFBlock> o) {
  assert(block_ == nullptr);
  block_ = o;
}


unique_ptr<double[]> ParallelDF::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  // first thing is to find the node
  const int inode = get_node(i); 

  // ask for the data to inode
  if (inode == mpi__->rank()) {
    return block_->get_block(i, id, j, jd, k, kd);
  } else {
    throw logic_error("ParallelDF::get_block is an intra-node function (or bug?)");
  }
  return unique_ptr<double[]>();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DFDist::add_direct_product(const vector<const double*> cd, const vector<const double*> dd, const double a) {
  if (cd.size() != dd.size()) throw logic_error("Illegal call of DFDist::DFDist");

  for (auto c = cd.begin(), d = dd.begin(); c != cd.end(); ++c, ++d)
    block_->add_direct_product(*c+block_->astart(), *d, a);
}


tuple<int, vector<shared_ptr<const Shell> > > DFDist::get_ashell(const vector<shared_ptr<const Shell> >& all) const {
  const int batchsize = naux_ / mpi__->size();
  const int start = batchsize*mpi__->rank();
  const int end   = mpi__->rank() != mpi__->size()-1 ? batchsize*mpi__->rank() + batchsize : naux_;

  int out1;
  vector<shared_ptr<const Shell> > out2;

  int num = 0;
  for (auto iter = all.begin(); iter != all.end(); ++iter) {
    if (num >= start && num < end) {
      if (out2.empty()) out1 = num;
      out2.push_back(*iter);
    }
    num += (*iter)->nbasis();
  }

  return tie(out1, out2);
}


void DFDist::common_init(const vector<shared_ptr<const Atom> >& atoms0, const vector<shared_ptr<const Atom> >& atoms1,
                         const vector<shared_ptr<const Atom> >& aux_atoms, const double throverlap, const bool compute_inverse) {

  Timer time;

  // 3index Integral is now made in DFBlock.
  vector<shared_ptr<const Shell> > ashell, b1shell, b2shell;
  for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms1) b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms0) b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

  // Decide how we distribute (dynamic distribution).
  // TODO we need a parallel queue server!

  // construction of DFBlock computes integrals
#ifndef HAVE_MPI_H
  block_ = shared_ptr<DFBlock>(new DFBlock(ashell, b1shell, b2shell, 0, 0, 0));
#else
  int astart;
  vector<shared_ptr<const Shell> > myashell;
  tie(astart, myashell) = get_ashell(ashell);
  block_ = shared_ptr<DFBlock>(new DFBlock(myashell, b1shell, b2shell, astart, 0, 0));
#endif

  time.tick_print("3-index ints");

  // make a global hash table
  make_table(mpi__->size());

  // generates a task of integral evaluations
  vector<DFIntTask_OLD<DFDist> > tasks;
  tasks.reserve(ashell.size()*ashell.size());
  data2_ = shared_ptr<ParaMatrix>(new ParaMatrix(naux_, naux_));

  int tmpa = 0;
  vector<int> aof;
  for (auto& i : ashell) { aof.push_back(tmpa); tmpa += i->nbasis(); }
  const shared_ptr<const Shell> b3(new Shell(atoms0.front()->shells().front()->spherical()));

  // naive static distribution
  int u = 0;
  auto o0 = aof.begin();
  for (auto& b0 : ashell) {
    auto o1 = aof.begin();
    for (auto& b1 : ashell) {
      if (*o0 <= *o1 && (u++ % mpi__->size() == mpi__->rank()))
        tasks.push_back(DFIntTask_OLD<DFDist>(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, vector<int>{*o0, *o1}, this));
      ++o1;
    }
    ++o0;
  }

  // these shell loops will be distributed across threads
  TaskQueue<DFIntTask_OLD<DFDist> > tq(tasks);
  tq.compute(resources__->max_num_threads());
  data2_->allreduce();

  time.tick_print("2-index ints");

  if (compute_inverse) {
    data2_->inverse_half(throverlap);
    time.tick_print("computing inverse");
  }
}


void DFDist::make_table(const int nblock) {
  unique_ptr<int[]> send(new int[2]);
  unique_ptr<int[]> rec(new int[2*mpi__->size()]);
  fill_n(rec.get(), 2*mpi__->size(), 0);

  send[0] = block_->astart();
  send[1] = block_->asize();
  mpi__->allgather(send.get(), 2, rec.get(), 2); 

  // reformatting to map
  int* buf = rec.get();
  for (int inode = 0; inode != mpi__->size(); ++inode, buf += 2) {
    global_table_.insert(make_pair(*buf+*(buf+1), inode));
    atable_.push_back(make_pair(*buf, *(buf+1)));
  }
}


pair<const double*, shared_ptr<RysInt> > DFDist::compute_batch(array<shared_ptr<const Shell>,4>& input) {
#ifdef LIBINT_INTERFACE
  shared_ptr<Libint> eribatch(new Libint(input));
#else
  shared_ptr<ERIBatch> eribatch(new ERIBatch(input, 2.0));
#endif
  eribatch->compute();
  return make_pair(eribatch->data(), eribatch);
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const Matrix> den) const {
  return compute_Jop(shared_from_this(), den);
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const ParallelDF> o, const shared_ptr<const Matrix> den) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  unique_ptr<double[]> tmp0 = o->compute_cd(den, data2_);
  // then compute J operator J_{rs} = |E*) (E|rs)
  shared_ptr<ParaMatrix> out = block_->form_mat(tmp0.get()+block_->astart());
  // all reduce
  out->allreduce();
  return out;
}


unique_ptr<double[]> ParallelDF::compute_cd(const shared_ptr<const Matrix> den, shared_ptr<const Matrix> dat2) const {
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  unique_ptr<double[]> tmp0(new double[naux_]);
  unique_ptr<double[]> tmp1(new double[naux_]);
  fill_n(tmp0.get(), naux_, 0.0);

  // D = (D|rs)*d_rs
  unique_ptr<double[]> tmp = block_->form_vec(den);
  copy_n(tmp.get(), block_->asize(), tmp0.get()+block_->astart());
  // All reduce
  mpi__->allreduce(tmp0.get(), naux_);

  // C = S^-1_CD D 
  dgemv_("N", naux_, naux_, 1.0, dat2->data(), naux_, tmp0.get(), 1, 0.0, tmp1.get(), 1);
  dgemv_("N", naux_, naux_, 1.0, dat2->data(), naux_, tmp1.get(), 1, 0.0, tmp0.get(), 1);
  return tmp0;
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform(const double* c, const size_t nocc) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(shared_from_this(), nocc));
  out->add_block(block_->transform_second(c, nocc));
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFHalfDist::compute_second_transform(const double* c, const size_t nocc) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nocc));
  out->add_block(block_->transform_third(c, nocc));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::copy() const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nindex1_));
  out->add_block(block_->copy());
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::clone() const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nindex1_));
  out->add_block(block_->clone());
  return out;
}


shared_ptr<DFDist> DFHalfDist::back_transform(const double* c) const{
  shared_ptr<DFDist> out(new DFDist(df_));
  out->add_block(block_->transform_second(c, df_->nindex1(), true));
  return out;
}


void DFHalfDist::rotate_occ(const double* d) {
  block_ = block_->transform_second(d, nindex1_);
}


shared_ptr<DFHalfDist> DFHalfDist::apply_density(const double* den) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nindex1_)); 
  out->add_block(block_->transform_third(den, nindex2_));
  return out;
}


shared_ptr<Matrix> DFHalfDist::compute_Kop_1occ(const double* den, const double a) const {
  return apply_density(den)->form_2index(df_, a);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFFullDist::copy() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->copy());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::clone() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->clone());
  return out;
}


void DFFullDist::symmetrize() {
  block_->symmetrize();
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DFHalfDist> DFFullDist::back_transform(const double* c) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nindex1_));
  out->add_block(block_->transform_third(c, df_->nindex2(), true));
  return out;
}


// 2RDM contractions
shared_ptr<DFFullDist> DFFullDist::apply_closed_2RDM() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->apply_rhf_2RDM());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->apply_uhf_2RDM(amat, bmat));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->apply_2RDM(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nindex1_, nindex2_));
  out->add_block(block_->apply_2RDM(rdm));
  return out;
}


shared_ptr<Matrix> DFFullDist::form_aux_2index_apply_J(const shared_ptr<const DFFullDist> o, const double a) const {
  shared_ptr<Matrix> tmp = ParallelDF::form_aux_2index(o, a);
  return shared_ptr<Matrix>(new Matrix(*tmp * *df_->data2()));
}


unique_ptr<double[]> DFFullDist::form_4index_1fixed(const shared_ptr<const DFFullDist> o, const double a, const size_t n) const {
  unique_ptr<double[]> out = block_->form_4index_1fixed(o->block_, a, n);
  const size_t size = block_->b1size() * block_->b2size() * o->block_->b1size();
  mpi__->allreduce(out.get(), size);
  return out;
}


void DFFullDist::set_product(const shared_ptr<const DFFullDist> o, const unique_ptr<double[]>& c, const int jdim, const size_t off) {
  block_->copy_block(o->block_->form_Dj(c, jdim), jdim, off*block_->asize());
}


//////// apply J functions ////////

shared_ptr<DFFullDist> DFFullDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFFullDist> out = clone();
#ifdef HAVE_MPI_H
  Timer mult(3);
  shared_ptr<DFDistT> work(new DFDistT(shared_from_this()));
  mult.tick_print("Form DFDistT");
  work = work->apply_J(d);
  mult.tick_print("Application of Inverse");
  work->get_paralleldf(out);
  mult.tick_print("Return DFDist");
#else
  out->block_->zero();
  out->block_->contrib_apply_J(block_, d);
#endif
  return out; 
}


shared_ptr<DFHalfDist> DFHalfDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFHalfDist> out = clone();
#ifdef HAVE_MPI_H
  Timer mult(3);
  shared_ptr<DFDistT> work(new DFDistT(shared_from_this()));
  mult.tick_print("Form DFDistT");
  work = work->apply_J(d);
  mult.tick_print("Applicatoin of Inverse");
  work->get_paralleldf(out);
  mult.tick_print("Return DFDist");
#else
  out->block_->zero();
  out->block_->contrib_apply_J(block_, d);
#endif
  return out;
}


