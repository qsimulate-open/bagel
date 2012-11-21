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
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <list>
#include <src/rysint/eribatch.h>
#include <src/rysint/libint.h>
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/df/df.h>

#include <src/df/dfinttask_old.h>

using namespace std;
using namespace chrono;
using namespace bagel;


ParallelDF::ParallelDF() {

}


unique_ptr<double[]> ParallelDF::form_2index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (blocks_.size() != o->blocks_.size()) throw logic_error("illegal call of ParallelDF::form_2index");
  const size_t size = blocks_.front()->b2size()*o->blocks_.front()->b2size();
  unique_ptr<double[]> out(new double[size]);
  fill_n(out.get(), size, 0.0);

  // loop over blocks
  for (auto i = blocks_.begin(), j = o->blocks_.begin(); i != blocks_.end(); ++i, ++j) { 
    unique_ptr<double[]> tmp = (!swap) ? (*i)->form_2index(*j, a) : (*j)->form_2index(*i, a);
    // accumulate
    daxpy_(size, 1.0, tmp, 1, out, 1);
  }
  return out;
}


unique_ptr<double[]> ParallelDF::form_4index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (blocks_.size() != o->blocks_.size()) throw logic_error("illegal call of ParallelDF::form_4index");
  const size_t size = blocks_.front()->b2size()*o->blocks_.front()->b2size() * blocks_.front()->b1size()*o->blocks_.front()->b1size();
  unique_ptr<double[]> out(new double[size]);
  fill_n(out.get(), size, 0.0);

  // loop over blocks
  for (auto i = blocks_.begin(), j = o->blocks_.begin(); i != blocks_.end(); ++i, ++j) { 
    unique_ptr<double[]> tmp = (!swap) ? (*i)->form_4index(*j, a) : (*j)->form_4index(*i, a);
    // accumulate
    daxpy_(size, 1.0, tmp, 1, out, 1);
  }
  return out;
}


shared_ptr<Matrix> ParallelDF::form_aux_2index(shared_ptr<const ParallelDF> o, const double a) const {
  // first allocate memory...
  const size_t idim = blocks_.back()->astart() + blocks_.back()->asize();
  const size_t jdim = o->blocks_.back()->astart() + o->blocks_.back()->asize();
  shared_ptr<Matrix> out(new Matrix(idim, jdim));

  // TODO to be distributed
  for (auto& j : o->blocks_)
    for (auto& i : blocks_)
      out->copy_block(i->astart(), j->astart(), i->asize(), j->asize(), i->form_aux_2index(j, a));
  return out;
}


void ParallelDF::daxpy(const double a, const shared_ptr<const ParallelDF> o) {
  if (blocks_.size() != o->blocks_.size()) throw logic_error("illegal call of ParallelDF::daxpy");
  auto ob = o->blocks_.begin();
  for (auto& i : blocks_) {
    i->daxpy(a, *ob);
    ++ob;
  }
}


void ParallelDF::scale(const double a) {
  for (auto& i : blocks_)
    i->scale(a);
}


void ParallelDF::add_block(shared_ptr<DFBlock> o) {
  blocks_.push_back(o);
}


void ParallelDF::make_table(const int inode) {
  // TODO need to broadcast
  for (auto& i : blocks_)
    global_table_.insert(make_pair(i->astart()+i->asize(), inode));
}


unique_ptr<double[]> ParallelDF::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  // first thing is to find the node
  const int inode = global_table_.upper_bound(i)->first;
  // now we still have only inode == 0 in the table
  const int mynode = 0;

  // ask for the data to inode
  if (inode == mynode) {
    for (auto& b : blocks_)
      if (b->astart() <= i && i < b->astart()+b->asize())
        return b->get_block(i, id, j, jd, k, kd);
    assert(false);
  } else {
    throw logic_error("not yet implemented ParallelDF::get_block");
  }
  return unique_ptr<double[]>();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DFDist::add_direct_product(const vector<const double*> cd, const vector<const double*> dd, const double a) {
  if (cd.size() != dd.size()) throw logic_error("Illegal call of DFDist::DFDist");
  for (auto& i : blocks_)
    for (auto c = cd.begin(), d = dd.begin(); c != cd.end(); ++c, ++d)
      i->add_direct_product(*c+i->astart(), *d, a);
}


void DFDist::common_init(const vector<shared_ptr<const Atom> >& atoms0, const vector<shared_ptr<const Atom> >& atoms1,
                         const vector<shared_ptr<const Atom> >& aux_atoms, const double throverlap, const bool compute_inverse) {

  auto tp0 = high_resolution_clock::now();

  // 3index Integral is now made in DFBlock.
  vector<shared_ptr<const Shell> > ashell, b1shell, b2shell;
  for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms1) b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms0) b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

  // Decide how we distribute (dynamic distribution).
  // TODO we need a parallel queue server!

  const int inode = 0;

  // construction of DFBlock computes integrals
#if 1
  blocks_.push_back(shared_ptr<DFBlock>(new DFBlock(ashell, b1shell, b2shell, 0, 0, 0)));
#else
  int cnt = 0;
  for (auto& i : ashell) {
    blocks_.push_back(shared_ptr<DFBlock>(new DFBlock(vector<shared_ptr<const Shell> >{i}, b1shell, b2shell, cnt, 0, 0)));
    cnt += i->nbasis();
  }
#endif

  // make a global hash table
  make_table(inode);

  // generates a task of integral evaluations
  vector<DFIntTask_OLD<DFDist> > tasks;
  data2_ = shared_ptr<Matrix>(new Matrix(naux_, naux_));

  int tmpa = 0;
  vector<int> aof;
  for (auto& i : ashell) { aof.push_back(tmpa); tmpa += i->nbasis(); }
  const shared_ptr<const Shell> b3(new Shell(atoms0.front()->shells().front()->spherical()));

  auto o0 = aof.begin();
  for (auto& b0 : ashell) {
    auto o1 = aof.begin();
    for (auto& b1 : ashell) {
      if (*o0 <= *o1)
        tasks.push_back(DFIntTask_OLD<DFDist>(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, vector<int>{*o0, *o1}, this));
      ++o1;
    }
    ++o0;
  }

  // these shell loops will be distributed across threads
  TaskQueue<DFIntTask_OLD<DFDist> > tq(tasks);
  tq.compute(resources__->max_num_threads());

  auto tp1 = high_resolution_clock::now();
  cout << "       - time spent for integral evaluation  " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl;

  if (compute_inverse) data2_->inverse_half(throverlap);
  auto tp2 = high_resolution_clock::now();
  cout << "       - time spent for computing inverse    " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp2-tp1).count()*0.001 << endl;

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


unique_ptr<double[]> DFDist::compute_Jop(const double* den) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  unique_ptr<double[]> tmp0 = compute_cd(den);
  unique_ptr<double[]> out(new double[nbasis0_*nbasis1_]);
  fill_n(out.get(), nbasis0_*nbasis1_, 0.0);
  // then compute J operator J_{rs} = |E*) (E|rs)
  for (auto& i : blocks_) {
    unique_ptr<double[]> tmp = i->form_mat(tmp0.get()+i->astart());
    daxpy_(nbasis0_*nbasis1_, 1.0, tmp, 1, out, 1);
  }
  return out;
}


unique_ptr<double[]> DFDist::compute_cd(const double* den) const {
  unique_ptr<double[]> tmp0(new double[naux_]);
  unique_ptr<double[]> tmp1(new double[naux_]);
  // D = (D|rs)*d_rs
  for (auto& i : blocks_) {
    unique_ptr<double[]> tmp = i->form_vec(den);
    copy_n(tmp.get(), i->asize(), tmp0.get()+i->astart());
  }
  // C = S^-1_CD D 
  dgemv_("N", naux_, naux_, 1.0, data2_->data(), naux_, tmp0.get(), 1, 0.0, tmp1.get(), 1);
  dgemv_("N", naux_, naux_, 1.0, data2_->data(), naux_, tmp1.get(), 1, 0.0, tmp0.get(), 1);
  return tmp0;
}


shared_ptr<DFHalfDist> DFDist::compute_half_transform(const double* c, const size_t nocc) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(shared_from_this(), nocc));
  for (auto& i : blocks_)
    out->add_block(i->transform_second(c, nocc));
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFHalfDist::compute_second_transform(const double* c, const size_t nocc) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc_, nocc));
  for (auto& i : blocks_)
    out->add_block(i->transform_third(c, nocc));
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::copy() const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nocc_));
  for (auto& i : blocks_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFHalfDist> DFHalfDist::clone() const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nocc_));
  for (auto& i : blocks_)
    out->add_block(i->clone());
  return out;
}


shared_ptr<DFDist> DFHalfDist::back_transform(const double* c) const{
  shared_ptr<DFDist> out(new DFDist(df_));
  for (auto& i : blocks_)
    out->add_block(i->transform_second(c, df_->nbasis1(), true));
  return out;
}


void DFHalfDist::rotate_occ(const double* d) {
  for (auto& i : blocks_)
    i = i->transform_second(d, nocc_);
}


shared_ptr<DFHalfDist> DFHalfDist::apply_density(const double* den) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nocc_)); 
  for (auto& i : blocks_)
    out->add_block(i->transform_third(den, nbasis_));
  return out;
}


unique_ptr<double[]> DFHalfDist::compute_Kop_1occ(const double* den) const {
  return apply_density(den)->form_2index(df_, 1.0);
}


shared_ptr<DFHalfDist> DFHalfDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFHalfDist> out = clone();
  for (auto& i : out->blocks_) {
    i->zero();
    for (auto& j : blocks_)
      i->contrib_apply_J(j, d);
  }
  return out; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


shared_ptr<DFFullDist> DFFullDist::copy() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->copy());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::clone() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->clone());
  return out;
}


void DFFullDist::symmetrize() {
  for (auto& i : blocks_)
    i->symmetrize();
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DFHalfDist> DFFullDist::back_transform(const double* c) const {
  shared_ptr<DFHalfDist> out(new DFHalfDist(df_, nocc1_));
  for (auto& i : blocks_)
    out->add_block(i->transform_third(c, df_->nbasis0(), true));
  return out;
}


// 2RDM contractions
shared_ptr<DFFullDist> DFFullDist::apply_closed_2RDM() const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->apply_rhf_2RDM());
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->apply_uhf_2RDM(amat, bmat));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->apply_2RDM(rdm, rdm1, nclosed, nact));
  return out;
}


shared_ptr<DFFullDist> DFFullDist::apply_2rdm(const double* rdm) const {
  shared_ptr<DFFullDist> out(new DFFullDist(df_, nocc1_, nocc2_));
  for (auto& i : blocks_)
    out->add_block(i->apply_2RDM(rdm));
  return out;
}


shared_ptr<Matrix> DFFullDist::form_aux_2index_apply_J(const shared_ptr<const DFFullDist> o, const double a) const {
  shared_ptr<Matrix> tmp = ParallelDF::form_aux_2index(o, a);
  return shared_ptr<Matrix>(new Matrix(*tmp * *df_->data2_));
}


unique_ptr<double[]> DFFullDist::form_4index_1fixed(const shared_ptr<const DFFullDist> o, const double a, const size_t n) const {
  const size_t size = blocks_.front()->b1size() * blocks_.front()->b2size() * o->blocks_.front()->b1size();
  unique_ptr<double[]> out(new double[size]);
  fill_n(out.get(), size, 0.0);

  // TODO will be distributed
  for (auto i = blocks_.begin(), j = o->blocks_.begin(); i != blocks_.end(); ++i, ++j) { 
    unique_ptr<double[]> tmp = (*i)->form_4index_1fixed(*j, a, n);
    daxpy_(size, 1.0, tmp, 1, out, 1);
  }
  return out;
}


void DFFullDist::set_product(const shared_ptr<const DFFullDist> o, const unique_ptr<double[]>& c, const int jdim, const size_t off) {
  auto j = o->blocks_.begin();
  for (auto& i : blocks_) {
    i->copy_block((*j)->form_Dj(c, jdim), jdim, off*i->asize());
    ++j;
  }
}


shared_ptr<DFFullDist> DFFullDist::apply_J(const shared_ptr<const Matrix> d) const {
  shared_ptr<DFFullDist> out = clone();
  for (auto& i : out->blocks_) {
    i->zero();
    for (auto& j : blocks_)
      i->contrib_apply_J(j, d);
  }
  return out; 
}

