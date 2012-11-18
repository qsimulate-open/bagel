//
// BAGEL - Parallel electron correlation program.
// Filename: dfdist.cc
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
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/df/dfdist.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <list>

using namespace std;
using namespace chrono;
using namespace bagel;


void DFDist::common_init(const vector<shared_ptr<const Atom> >& atoms0, const vector<shared_ptr<const Atom> >& atoms1,
                         const vector<shared_ptr<const Atom> >& aux_atoms, const double throverlap, const bool compute_inverse) {

#if 0
  // this will be distributed in the future.
  auto tp0 = high_resolution_clock::now();

  // 3index Integral is now made in DFBlock.
  vector<shared_ptr<const Shell> > ashell, b1shell, b2shell;
  for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms1) b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms0) b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

  // Decide how we distribute (dynamic distribution).
  // TODO we need a parallel queue server!
  data_ = shared_ptr<DFBlock>(new DFBlock(ashell, b1shell, b2shell, 0, 0, 0));

  // generates a task of integral evaluations
  vector<DFIntTask_OLD> tasks;
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
        tasks.push_back(DFIntTask_OLD(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, vector<int>{*o0, *o1}, this));
      ++o1;
    }
    ++o0;
  }

  // these shell loops will be distributed across threads
  TaskQueue<DFIntTask_OLD> tq(tasks);
  tq.compute(resources__->max_num_threads());
  auto tp1 = high_resolution_clock::now();
  cout << "       - time spent for integral evaluation  " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl;

  if (compute_inverse) data2_->inverse_half(throverlap);
  auto tp2 = high_resolution_clock::now();
  cout << "       - time spent for computing inverse    " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp2-tp1).count()*0.001 << endl;
#endif

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
  fill(out.get(), out.get()+nbasis0_*nbasis1_, 0.0);
  // then compute J operator J_{rs} = |E*) (E|rs)
  for (auto& i : blocks_)
    dgemv_("T", i->asize(), nbasis0_*nbasis1_, 1.0, data_->get(), i->asize(), tmp0.get()+i->astart(), 1, 1.0, out.get(), 1);
  return out;
}


unique_ptr<double[]> DFDist::compute_cd(const double* den) const {
  unique_ptr<double[]> tmp0(new double[naux_]);
  unique_ptr<double[]> tmp1(new double[naux_]);
  // D = (D|rs)*d_rs
  for (auto& i : blocks_)
    dgemv_("N", i->asize(), nbasis0_*nbasis1_, 1.0, data_->get(), i->asize(), den, 1, 0.0, tmp0.get()+i->astart(), 1);
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


void DFFullDist::daxpy(const double a, const DFFullDist& o) {
  if (blocks_.size() != o.blocks_.size()) throw logic_error("illegal call of DFFullDist::daxpy");
  auto ob = o.blocks_.begin();
  for (auto& i : blocks_) {
    i->daxpy(a, *ob);
    ++ob;
  }
}


void DFFullDist::daxpy(const double a, const DFHalfDist& o) {
  if (blocks_.size() != o.blocks_.size()) throw logic_error("illegal call of DFFullDist::daxpy");
  auto ob = o.blocks_.begin();
  for (auto& i : blocks_) {
    i->daxpy(a, *ob);
    ++ob;
  }
}


void DFFullDist::scale(const double a) {
  for (auto& i : blocks_)
    i->scale(a);
}


void DFFullDist::symmetrize() {
  for (auto& i : blocks_)
    i->symmetrize();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
DF_AO::DF_AO(const int nbas0, const int nbas1, const int naux, const vector<const double*> cd, const vector<const double*> dd)
 : DensityFit(nbas0, nbas1, naux) {
  assert(cd.size() == dd.size());

  // initialize to zero
  unique_ptr<double[]> buf(new double[naux_*nbasis0_*nbasis1_]);
  fill(buf.get(), buf.get()+size(), 0.0);

  for (auto citer = cd.begin(), diter = dd.begin(); citer != cd.end(); ++citer, ++diter) {
    dger_(naux_, nbasis0_*nbasis1_, 1.0, *citer, 1, *diter, 1, buf.get(), naux_);
  }
  data_ = shared_ptr<DFBlock>(new DFBlock(buf, naux, nbas1, nbas0, 0,0,0));
}


shared_ptr<DF_Half> DF_Half::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
  dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, d->data2_->data(), naux_, data_->get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, out));
}


shared_ptr<DF_Half> DF_Half::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  unique_ptr<double[]> jj(new double[naux_*naux_]);
  dgemm_("N", "N", naux_, naux_, naux_, 1.0, d->data2_->data(), naux_, d->data2_->data(), naux_, 0.0, jj.get(), naux_);

  unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
  dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, jj.get(), naux_, data_->get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, out));
}


void DF_Half::form_2index(unique_ptr<double[]>& target, const double a, const double b) const {
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, nbasis_, common, a, data_->get(), common, data_->get(), common, b, target.get(), nbasis_);
}


unique_ptr<double[]> DF_Half::form_2index(shared_ptr<const DF_Full> o, const double a, const double b) const {
  assert(b == 0);
  unique_ptr<double[]> tmp(new double[nbasis_*o->nocc2()]);
  form_2index(tmp, o, a, b);
  return tmp;
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DF_Full> o, const double a, const double b) const {
  if (nocc_ != o->nocc1()) throw logic_error("nocc_ and o->nocc1() should be the same: in DF_Half::form_2index");
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, o->nocc2(), common, a, data_->get(), common, o->data_->get(), common, b, target.get(), nbasis_);
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DensityFit> o, const double a) const {
  fill(target.get(), target.get()+nocc_*nbasis_, 0.0);
  for (size_t i = 0; i != nbasis_; ++i)
    dgemm_("T", "N", nocc_, nbasis_, naux_, a, data_->get()+i*naux_*nocc_, naux_, o->data_->get()+i*naux_*nbasis_, naux_,
                                            1.0, target.get(), nocc_);
}
unique_ptr<double[]> DF_Half::form_2index(shared_ptr<const DensityFit> o, const double a) const {
  unique_ptr<double[]> out(new double[nocc_*nbasis_]);
  form_2index(out, o, a);
  return move(out);
}


void DF_Half::form_4index(unique_ptr<double[]>& target) const {
  const int ndim = nbasis_ * nocc_;
  dgemm_("T", "N", ndim, ndim, naux_, 1.0, data_->get(), naux_, data_->get(), naux_, 0.0, target.get(), ndim);
}

unique_ptr<double[]> DF_Half::form_4index() const {
  const size_t ndim = nbasis_ * nocc_;
  unique_ptr<double[]> out(new double[ndim*ndim]);
  form_4index(out);
  return move(out);
}


unique_ptr<double[]> DF_Half::compute_Kop_1occ(const double* den) const {
  return apply_density(den)->form_2index(df_);
}


unique_ptr<double[]> DF_Half::form_aux_2index(const shared_ptr<const DF_Half> o) const {
  unique_ptr<double[]> out(new double[naux_*naux_]);
  if (nocc_*nbasis_ != o->nocc_*o->nbasis_) throw logic_error("wrong call to DF_Full::form_aux_2index");
  dgemm_("N", "T", naux_, naux_, nocc_*nbasis_, 1.0, data_->get(), naux_, o->data_->get(), naux_, 0.0, out.get(), naux_);
  return out;
}



shared_ptr<DF_Half> DF_Half::apply_density(const double* den) const {
  unique_ptr<double[]> buf(new double[naux_*nbasis_*nocc_]);
  dgemm_("N", "N", naux_*nocc_, nbasis_, nbasis_, 1.0, data_->get(), naux_*nocc_, den, nbasis_, 0.0, buf.get(), naux_*nocc_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, buf));
}


void DF_Half::rotate_occ(const double* d) {
  unique_ptr<double[]> buf(new double[naux_*nocc_]);
  for (int i = 0; i != nbasis_; ++i) {
    dgemm_("N", "N", naux_, nocc_, nocc_, 1.0, data_->get()+i*naux_*nocc_, naux_, d, nocc_, 0.0, buf.get(), naux_);
    copy_n(buf.get(), naux_*nocc_, data_->get()+i*naux_*nocc_); 
  }
}


shared_ptr<DF_Full> DF_Full::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, d->data2_->data(), naux_, data_->get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  unique_ptr<double[]> jj(new double[naux_*naux_]);
  dgemm_("N", "N", naux_, naux_, naux_, 1.0, d->data2_->data(), naux_, d->data2_->data(), naux_, 0.0, jj.get(), naux_);

  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, jj.get(), naux_, data_->get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm) const {
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "T", naux_, nocc1_*nocc2_, nocc1_*nocc2_, 1.0, data_->get(), naux_, rdm, nocc1_*nocc2_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == nocc1_ && nocc1_ == nocc2_);
  // checking if natural orbitals...
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__) throw logic_error("DF_Full::apply_2rdm should be called with natural orbitals");
  }
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  fill(out.get(), out.get()+nocc1_*nocc2_*naux_, 0.0);
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      daxpy_(naux_, -2.0, data_->get()+naux_*(j+nocc1_*i), 1, out.get()+naux_*(j+nocc1_*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 1.0, data_->get()+naux_*(i+nocc1_*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 4.0, diagsum.get(), 1, out.get()+naux_*(i+nocc1_*i), 1);

  // act-act part
  // compress
  unique_ptr<double[]> buf(new double[nact*nact*naux_]);
  unique_ptr<double[]> buf2(new double[nact*nact*naux_]);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(naux_, data_->get()+naux_*(j+nclosed+nocc1_*(i+nclosed)), 1, buf.get()+naux_*(j+nact*i),1);
  // multiply
  dgemm_("N", "N", naux_, nact*nact, nact*nact, 1.0, buf.get(), naux_, rdm, nact*nact, 0.0, buf2.get(), naux_);
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(naux_, buf2.get()+naux_*(j+nact*i),1, out.get()+naux_*(j+nclosed+nocc1_*(i+nclosed)), 1);

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  for (int i = 0; i != nact; ++i)
    daxpy_(naux_, 2.0*rdm1[i+nact*i], diagsum.get(), 1, out.get()+naux_*(i+nclosed+nocc1_*(i+nclosed)), 1);
  unique_ptr<double[]> diagsum2(new double[naux_]);
  dgemv_("N", naux_, nact*nact, 1.0, buf.get(), naux_, rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 2.0, diagsum2.get(), 1, out.get()+naux_*(i+nocc1_*i), 1);
  // exchange contribution
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nclosed; ++j) {
      daxpy_(naux_, -rdm1[i+nact*i], data_->get()+naux_*(j+nocc1_*(i+nclosed)), 1, out.get()+naux_*(j+nocc1_*(i+nclosed)), 1);
      daxpy_(naux_, -rdm1[i+nact*i], data_->get()+naux_*(i+nclosed+nocc1_*j), 1, out.get()+naux_*(i+nclosed+nocc1_*j), 1);
    }
  }

  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


// forms all-internal 4-index MO integrals
void DF_Full::form_4index(unique_ptr<double[]>& target) const {
  const int dim = nocc1_ * nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, dim, naux, 1.0, data_->get(), naux, data_->get(), naux, 0.0, target.get(), dim);
}


unique_ptr<double[]> DF_Full::form_4index() const {
  unique_ptr<double[]> out(new double[nocc1_ * nocc2_ * nocc1_ * nocc2_]);
  form_4index(out);
  return move(out);
}


void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nocc1_ * o->nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, odim, naux, 1.0, data_->get(), naux, o->data_->get(), naux, 0.0, target.get(), dim);
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DF_Full> o) const {
  unique_ptr<double[]> out(new double[nocc1_ * nocc2_ * o->nocc1_ * o->nocc2_]);
  form_4index(out, o);
  return move(out);
}


// for MP2-like quantities
void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o, const size_t n) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nocc1_; // o->nocc2_ is fixed at n;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, odim, naux, 1.0, data_->get(), naux, o->data_->get()+naux*odim*n, naux, 0.0, target.get(), dim);
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DF_Full> o, const size_t n) const {
  unique_ptr<double[]> target(new double[o->nocc1_*nocc1_*nocc2_]);
  form_4index(target, o, n);
  return move(target);
}


// Joperator. Note that (r,s) runs first; i.e., in the operator form
void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DensityFit> o) const {
  const size_t dim = nocc1_ * nocc2_;
  const size_t odim = o->nbasis0() * o->nbasis1();
  shared_ptr<DF_Full> tmp = this->apply_J();
  dgemm_("T", "N", odim, dim, naux_, 1.0, o->data_->get(), naux_, tmp->data_->get(), naux_, 0.0, target.get(), odim);
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DensityFit> o) const {
  const size_t dim0 = nocc1_ * o->nbasis0();
  const size_t dim1 = nocc2_ * o->nbasis1();
  assert(nocc1_ == nocc2_); // TODO <- do you need this?
  unique_ptr<double[]> out(new double[dim0*dim1]);
  form_4index(out, o);
  return move(out);
}


unique_ptr<double[]> DF_Full::form_aux_2index(const shared_ptr<const DF_Full> o) const {
  unique_ptr<double[]> out(new double[naux_*naux_]);
  if (nocc1_*nocc2_ != o->nocc1_*o->nocc2_) throw logic_error("wrong call to DF_Full::form_aux_2index");
  dgemm_("N", "T", naux_, naux_, nocc1_*nocc2_, 1.0, data_->get(), naux_, o->data_->get(), naux_, 0.0, out.get(), naux_);
  return out;
}


unique_ptr<double[]> DF_Full::form_aux_2index_apply_J(const shared_ptr<const DF_Full> o) const {
  unique_ptr<double[]> tmp(new double[naux_*naux_]);
  if (nocc1_*nocc2_ != o->nocc1_*o->nocc2_) throw logic_error("wrong call to DF_Full::form_aux_2index");
  dgemm_("N", "T", naux_, naux_, nocc1_*nocc2_, 1.0, data_->get(), naux_, o->data_->get(), naux_, 0.0, tmp.get(), naux_);
  unique_ptr<double[]> out(new double[naux_*naux_]);
  dgemm_("N", "N", naux_, naux_, naux_, 1.0, tmp.get(), naux_, df_->data2_->data(), naux_, 0.0, out.get(), naux_);
  return out;
}


void DF_Full::form_2index(unique_ptr<double[]>& target, const shared_ptr<const DF_Half> o, const double a) {
  assert(nocc1() == o->nocc());
  dgemm_("T", "N", nocc2_, o->nbasis(), nocc1_*naux_, a, data_->get(), nocc1_*naux_, o->data_->get(), nocc1_*naux_, 0.0, target.get(), nocc2_);
}


unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Half> o, const double a) {
  unique_ptr<double[]> out(new double[nocc2_*o->nbasis()]);
  form_2index(out, o, a);
  return out;
}


void DF_Full::form_2index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o, const double a) {
  assert(nocc1_ == o->nocc1_);
  dgemm_("T", "N", nocc2_, o->nocc2(), nocc1_*naux_, a, data_->get(), nocc1_*naux_, o->data_->get(), nocc1_*naux_, 0.0, target.get(), nocc2_);
}


unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Full> o, const double a) {
  unique_ptr<double[]> out(new double[nocc2_*o->nocc2()]);
  form_2index(out, o, a);
  return out;
}


void DF_Full::set_product(const shared_ptr<const DF_Full> o, const unique_ptr<double[]>& c, const int jdim, const size_t off) {
  dgemm_("N", "N", naux(), jdim, nocc1()*nocc2(), 1.0, o->data_->get(), naux(), c.get(), nocc1()*nocc2(), 0.0, data_->get()+off, naux());
}



shared_ptr<DF_Full> DF_Full::apply_closed_2RDM() const {
  assert(nocc1_ == nocc2_);
  const int nocc = nocc1_;
  unique_ptr<double[]> d(new double[size()]);
  fill(d.get(), d.get()+size(), 0.0);
  // exchange contributions
  daxpy_(size(), -2.0, data_->get(), 1, d.get(), 1);
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, 1.0, data_->get()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
  }
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, 4.0, diagsum.get(), 1, d.get()+naux_*(i+nocc*i), 1);
  }
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DF_Full> DF_Full::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(nocc1_ == nocc2_);
  const int nocc = nocc1_;
  unique_ptr<double[]> d(new double[size()]);
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data_->get(), naux_*nocc, amat, nocc, 0.0, d2.get(), naux_*nocc);
    for (int i = 0; i != nocc; ++i) {
      dgemm_("N", "N", naux_, nocc, nocc, -1.0, d2.get()+naux_*nocc*i, naux_, amat, nocc, 0.0, d.get()+naux_*nocc*i, naux_);
    }
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data_->get(), naux_*nocc, bmat, nocc, 0.0, d2.get(), naux_*nocc);
    for (int i = 0; i != nocc; ++i) {
      dgemm_("N", "N", naux_, nocc, nocc, -1.0, d2.get()+naux_*nocc*i, naux_, bmat, nocc, 1.0, d.get()+naux_*nocc*i, naux_);
    }
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc];
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, sum[i], data_->get()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
  }
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, sum[i], diagsum.get(), 1, d.get()+naux_*(i+nocc*i), 1);
  }
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DF_Half> DF_Full::back_transform(const double* c) const{
  const int nbas = df_->nbasis1();
  unique_ptr<double[]> d(new double[nbas*nocc1_*naux_]);
  dgemm_("N", "T", naux_*nocc1_, nbas, nocc2_, 1.0, data_->get(), naux_*nocc1_, c, nbas, 0.0, d.get(), naux_*nocc1_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc1_, d));
}


shared_ptr<DF_AO> DF_Half::back_transform(const double* c) const{
  const int nbas = df_->nbasis0();
  unique_ptr<double[]> d(new double[nbas*nbasis_*naux_]);
  for (int i = 0; i != nbasis_; ++i)
    dgemm_("N", "T", naux_, nbas, nocc_, 1.0, data_->get()+i*naux_*nocc_, naux_, c, nbas, 0.0, d.get()+i*naux_*nbas, naux_);
  return shared_ptr<DF_AO>(new DF_AO(nbas, nbasis_, naux_, d));
}
#endif


