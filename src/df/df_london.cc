//
// BAGEL - Parallel electron correlation program.
// Filename: df_london.cc
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


#include <src/df/df_london.h>
#include <src/df/dfdistt_london.h>
#include <src/df/paralleldf_london.h>
#include <src/integral/rys/eribatch.h>

using namespace std;
using namespace bagel;
using namespace btas;


shared_ptr<DFDist_London> DFDist_London::copy() const {
  auto out = make_shared<DFDist_London>(df_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFDist_London> DFDist_London::clone() const {
  auto out = make_shared<DFDist_London>(df_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


void DFDist_London::add_direct_product(const vector<shared_ptr<const ZMatrix>> cd, const vector<shared_ptr<const ZMatrix>> dd, const double a) {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  if (cd.size() != dd.size()) throw logic_error("Illegal call of DFDist_London::DFDist_London");

  auto d = dd.begin();
  for (auto& c : cd) {
    shared_ptr<const ZMatrix> aslice = c->get_submatrix(block_[0]->astart(), 0, block_[0]->asize(), 1);
    block_[0]->add_direct_product(aslice, *d++, a);
  }
  assert(d == dd.end());
}


tuple<int, vector<shared_ptr<const Shell>>> DFDist_London::get_ashell(const vector<shared_ptr<const Shell>>& all) {
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


void DFDist_London::compute_2index(const vector<shared_ptr<const Shell>>& ashell, const double throverlap, const bool compute_inverse) {
  Timer time;

  // generates a task of integral evaluations
  TaskQueue<DFIntTask_OLD<DFDist_London>> tasks(ashell.size()*ashell.size());

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


shared_ptr<const StaticDist> DFDist_London::make_table(const size_t astart) {
  vector<size_t> rec(mpi__->size());
  fill(rec.begin(), rec.end(), 0);

  mpi__->allgather(&astart, 1, rec.data(), 1);
  rec.push_back(naux_);

  return make_shared<const StaticDist>(rec);
}


pair<const double*, shared_ptr<RysIntegral<double, Int_t::Standard>>> DFDist_London::compute_batch(array<shared_ptr<const Shell>,4>& input) {
  shared_ptr<RysIntegral<double, Int_t::Standard>> londoneribatch = make_shared<ERIBatch>(input, 2.0);
  londoneribatch->compute();
  return make_pair(londoneribatch->data(), londoneribatch);
}


shared_ptr<DFHalfDist_London> DFDist_London::compute_half_transform(shared_ptr<const View2<complex<double>>> c) const {
  const int nocc = c->range(1).size();
  auto out = make_shared<DFHalfDist_London>(shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_second(c));
  return out;
}


shared_ptr<DFHalfDist_London> DFDist_London::compute_half_transform_swap(shared_ptr<const View2<complex<double>>> c) const {
  const int nocc = c->range(1).size();
  auto out = make_shared<DFHalfDist_London>(shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c)->swap());
  return out;
}


// TODO will be deprecated
shared_ptr<DFHalfDist_London> DFDist_London::compute_half_transform(const shared_ptr<const ZMatrix> c) const {
  const int nocc = c->mdim();
  auto out = make_shared<DFHalfDist_London>(shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_second(c));
  return out;
}


// TODO will be deprecated
shared_ptr<DFHalfDist_London> DFDist_London::compute_half_transform_swap(const shared_ptr<const ZMatrix> c) const {
  const int nocc = c->mdim();
  auto out = make_shared<DFHalfDist_London>(shared_from_this(), nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c)->swap());
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist_London> DFHalfDist_London::compute_second_transform(shared_ptr<const View2<complex<double>>> c) const {
  const int nocc = c->range(1).size();
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c));
  return out;
}


// TODO will be deprecated
shared_ptr<DFFullDist_London> DFHalfDist_London::compute_second_transform(const shared_ptr<const ZMatrix> c) const {
  const int nocc = c->mdim();
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nocc);
  for (auto& i : block_)
    out->add_block(i->transform_third(c));
  return out;
}


shared_ptr<DFHalfDist_London> DFHalfDist_London::copy() const {
  auto out = make_shared<DFHalfDist_London> (df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFHalfDist_London> DFHalfDist_London::clone() const {
  auto out = make_shared<DFHalfDist_London>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


shared_ptr<DFDist_London> DFHalfDist_London::back_transform(shared_ptr<const View2<complex<double>>> c) const{
  assert(df_->nindex1() == c->range(0).size());
  auto out = make_shared<DFDist_London>(df_);
  for (auto& i : block_)
    out->add_block(i->transform_second(c, true));
  return out;
}


// TODO will be deprecated
shared_ptr<DFDist_London> DFHalfDist_London::back_transform(const shared_ptr<const ZMatrix> c) const{
  assert(df_->nindex1() == c->ndim());
  auto out = make_shared<DFDist_London>(df_);
  for (auto& i : block_)
    out->add_block(i->transform_second(c, true));
  return out;
}


void DFHalfDist_London::rotate_occ(const shared_ptr<const ZMatrix> d) {
  assert(nindex1_ == d->mdim());
  for (auto& i : block_)
    i = i->transform_second(d);
}


shared_ptr<DFHalfDist_London> DFHalfDist_London::apply_density(const shared_ptr<const ZMatrix> den) const {
  assert(den->mdim() == nindex2_);
  auto out = make_shared<DFHalfDist_London>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->transform_third(den));
  return out;
}


shared_ptr<ZMatrix> DFHalfDist_London::compute_Kop_1occ(const shared_ptr<const ZMatrix> den, const double a) const {
  return apply_density(den)->form_2index(df_, a);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist_London> DFFullDist_London::copy() const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFFullDist_London> DFFullDist_London::clone() const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->clone());
  return out;
}


void DFFullDist_London::rotate_occ1(const shared_ptr<const ZMatrix> d) {
  assert(nindex1_ == d->mdim());
  for (auto& i : block_)
    i = i->transform_second(d);
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DFHalfDist_London> DFFullDist_London::back_transform(shared_ptr<const View2<complex<double>>> c) const {
  assert(c->range(0).size() == df_->nindex2());
  auto out = make_shared<DFHalfDist_London>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->transform_third(c, true));
  return out;
}


// TODO will be deprecated
shared_ptr<DFHalfDist_London> DFFullDist_London::back_transform(const shared_ptr<const ZMatrix> c) const {
  assert(c->ndim() == df_->nindex2());
  auto out = make_shared<DFHalfDist_London>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->transform_third(c, true));
  return out;
}

/*
// 2RDM contractions
shared_ptr<DFFullDist_London> DFFullDist_London::apply_closed_2RDM(const double scale_exch) const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_rhf_2RDM(scale_exch));
  return out;
}


shared_ptr<DFFullDist_London> DFFullDist_London::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_uhf_2RDM(amat, bmat));
  return out;
}


shared_ptr<DFFullDist_London> DFFullDist_London::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_2RDM(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<DFFullDist_London> DFFullDist_London::apply_2rdm(const double* rdm) const {
  auto out = make_shared<DFFullDist_London>(df_, nindex1_, nindex2_);
  for (auto& i : block_)
    out->add_block(i->apply_2RDM(rdm));
  return out;
}
*/

shared_ptr<ZMatrix> DFFullDist_London::form_aux_2index_apply_J(const shared_ptr<const DFFullDist_London> o, const double a) const {
  shared_ptr<ZMatrix> tmp = ParallelDF_London::form_aux_2index(o, a);
  return make_shared<ZMatrix>(*tmp * *df_->data2());
}


shared_ptr<ZMatrix> DFFullDist_London::form_4index_1fixed(const shared_ptr<const DFFullDist_London> o, const double a, const size_t n) const {
  // TODO needs more work
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<ZMatrix> out = block_[0]->form_4index_1fixed(o->block_[0], a, n);
  if (!serial_)
    out->allreduce();
  return out;
}


void DFFullDist_London::add_product(const shared_ptr<const DFFullDist_London> o, const shared_ptr<const ZMatrix> c, const int jdim, const size_t off, const double fac) {
  // TODO needs more work
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  block_[0]->add_block(o->block_[0]->form_Dj(c, jdim), jdim, off*block_[0]->asize(), fac);
}


shared_ptr<DFFullDist_London> DFFullDist_London::swap() const {
  auto out = make_shared<DFFullDist_London>(shared_from_this(), nocc2(), nocc1());
  for (auto& i : block_)
    out->add_block(i->swap());
  return out;
}


//////// apply J functions ////////

shared_ptr<DFFullDist_London> DFFullDist_London::apply_J(const shared_ptr<const ZMatrix> d) const {
  shared_ptr<DFFullDist_London> out = clone();
#ifdef HAVE_MPI_H
  if (!serial_) {
    Timer mult(3);
    auto work = make_shared<DFDistT_London>(shared_from_this());
    mult.tick_print("Form DFDistT_London");
    work = work->apply_J(d);
    mult.tick_print("Application of Inverse");
    work->get_paralleldf(out);
    mult.tick_print("Return DFDist_London");
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


shared_ptr<DFHalfDist_London> DFHalfDist_London::apply_J(const shared_ptr<const ZMatrix> d) const {
  shared_ptr<DFHalfDist_London> out = clone();
#ifdef HAVE_MPI_H
  if (!serial_) {
    Timer mult(3);
    auto work = make_shared<DFDistT_London>(shared_from_this());
    mult.tick_print("Form DFDistT_London");
    work = work->apply_J(d);
    mult.tick_print("Application of Inverse");
    work->get_paralleldf(out);
    mult.tick_print("Return DFDist_London");
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


