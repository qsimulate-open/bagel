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
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <list>

#include <src/df/dfinttask_old.h>

using namespace std;
using namespace chrono;
using namespace bagel;

#ifndef USE_DFDIST


void DensityFit::common_init(const vector<shared_ptr<const Atom> >& atoms0, const vector<shared_ptr<const Atom> >& atoms1,
                             const vector<shared_ptr<const Atom> >& aux_atoms, const double throverlap, const bool compute_inverse) {

  // this will be distributed in the future.
  auto tp0 = high_resolution_clock::now();

  // 3index Integral is now made in DFBlock.
  vector<shared_ptr<const Shell> > ashell, b1shell, b2shell;
  for (auto& i : aux_atoms) ashell.insert(ashell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms1) b1shell.insert(b1shell.end(), i->shells().begin(), i->shells().end());
  for (auto& i : atoms0) b2shell.insert(b2shell.end(), i->shells().begin(), i->shells().end());

  // TODO in the future DFBlock should be devided into pieces that are to be distributed
  //      ... DFBlock only takes care of intra-node parallelization
  data_ = shared_ptr<DFBlock>(new DFBlock(ashell, b1shell, b2shell, 0, 0, 0));

  // generates a task of integral evaluations
  vector<DFIntTask_OLD<DensityFit> > tasks;
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
        tasks.push_back(DFIntTask_OLD<DensityFit>(array<shared_ptr<const Shell>,4>{{b1, b3, b0, b3}}, vector<int>{*o0, *o1}, this));
      ++o1;
    }
    ++o0;
  }

  // these shell loops will be distributed across threads
  TaskQueue<DFIntTask_OLD<DensityFit> > tq(tasks);
  tq.compute(resources__->max_num_threads());
  auto tp1 = high_resolution_clock::now();
  cout << "       - time spent for integral evaluation  " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl;

  if (compute_inverse) data2_->inverse_half(throverlap);
  auto tp2 = high_resolution_clock::now();
  cout << "       - time spent for computing inverse    " << setprecision(2) << setw(10) << duration_cast<milliseconds>(tp2-tp1).count()*0.001 << endl;

}


unique_ptr<double[]> DensityFit::compute_Jop(const double* den) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  unique_ptr<double[]> tmp = compute_cd(den);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return data_->form_mat(tmp.get());
}


unique_ptr<double[]> DensityFit::compute_cd(const double* den) const {
  unique_ptr<double[]> tmp0 = data_->form_vec(den);
  unique_ptr<double[]> tmp1(new double[naux_]);
  dgemv_("N", naux_, naux_, 1.0, data2_->data(), naux_, tmp0.get(), 1, 0.0, tmp1.get(), 1);
  dgemv_("N", naux_, naux_, 1.0, data2_->data(), naux_, tmp1.get(), 1, 0.0, tmp0.get(), 1);
  return tmp0;
}


shared_ptr<DF_Half> DensityFit::compute_half_transform(const double* c, const size_t nocc) const {
  return shared_ptr<DF_Half>(new DF_Half(shared_from_this(), nocc, data_->transform_second(c, nocc)));
}


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


shared_ptr<DF_Half> DF_Half::clone() const {
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, data_->clone()));
}


shared_ptr<DF_Half> DF_Half::copy() const {
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, data_->copy()));
}


shared_ptr<DF_Half> DF_Half::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  shared_ptr<DF_Half> out = clone();
  out->data_->zero();
  out->data_->contrib_apply_J(data_, d->data2_);
  return out;
}


shared_ptr<DF_Half> DF_Half::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  shared_ptr<DF_Half> out = clone();
  out->data_->zero();
  shared_ptr<const Matrix> jj(new Matrix(*d->data2_ * *d->data2_));
  out->data_->contrib_apply_J(data_, jj);
  return out;
}


shared_ptr<DF_Full> DF_Half::compute_second_transform(const double* c, const size_t nocc) const {
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc_, nocc, data_->transform_third(c, nocc)));
}


unique_ptr<double[]> DF_Half::compute_Kop_1occ(const double* den) const {
  return apply_density(den)->form_2index(df_, 1.0);
}



shared_ptr<DF_Half> DF_Half::apply_density(const double* den) const {
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, data_->transform_third(den, nbasis_)));
}


void DF_Half::rotate_occ(const double* d) { data_ = data_->transform_second(d, nocc_); }


shared_ptr<DF_Full> DF_Full::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  shared_ptr<DF_Full> out = clone();
  out->data_->zero();
  out->data_->contrib_apply_J(data_, d->data2_);
  return out;
}


shared_ptr<DF_Full> DF_Full::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->has_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  shared_ptr<DF_Full> out = clone();
  out->data_->zero();
  shared_ptr<const Matrix> jj(new Matrix(*d->data2_ * *d->data2_));
  out->data_->contrib_apply_J(data_, jj);
  return out;
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm) const {
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->apply_2RDM(rdm)));
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->apply_2RDM(rdm, rdm1, nclosed, nact)));
}



unique_ptr<double[]> DF_Half::form_4index(const shared_ptr<const DF_Half> o, const double a) const { return data_->form_4index(o->data_, a); }
unique_ptr<double[]> DF_Full::form_4index(const double a) const { return data_->form_4index(data_, a); }
unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DF_Full> o, const double a) const { return data_->form_4index(o->data_, a); }
// !! CAUTION !!
// Joperator. Note that (r,s) runs first; i.e., in the operator form
unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DensityFit> o, const double a) const { return o->data_->form_4index(data_, a); };


// for MP2-like quantities
unique_ptr<double[]> DF_Full::form_4index_1fixed(const shared_ptr<const DF_Full> o, const double a, const size_t n) const { return data_->form_4index_1fixed(o->data_, a, n); }


shared_ptr<Matrix> DF_Half::form_aux_2index(const shared_ptr<const DF_Half> o, const double a) const { return data_->form_aux_2index(o->data_, a); }
shared_ptr<Matrix> DF_Full::form_aux_2index(const shared_ptr<const DF_Full> o, const double a) const { return data_->form_aux_2index(o->data_, a); }

shared_ptr<Matrix> DF_Full::form_aux_2index_apply_J(const shared_ptr<const DF_Full> o, const double a) const {
  shared_ptr<Matrix> tmp = data_->form_aux_2index(o->data_, a);
  return shared_ptr<Matrix>(new Matrix(*tmp * *df_->data2_));
}



unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Half> o, const double a) { return data_->form_2index(o->data_, a); }
unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Full> o, const double a) { return data_->form_2index(o->data_, a); }
unique_ptr<double[]> DF_Half::form_2index(const shared_ptr<const DF_Half> o, const double a) const { return data_->form_2index(o->data_, a); }
unique_ptr<double[]> DF_Half::form_2index(const shared_ptr<const DF_Full> o, const double a) const { return data_->form_2index(o->data_, a); }
unique_ptr<double[]> DF_Half::form_2index(const shared_ptr<const DensityFit> o, const double a) const { return data_->form_2index(o->data_, a); }


void DF_Full::set_product(const shared_ptr<const DF_Full> o, const unique_ptr<double[]>& c, const int jdim, const size_t off) {
  data_->copy_block(o->data_->form_Dj(c, jdim), jdim, off);
}


shared_ptr<DF_Full> DF_Full::clone() const { return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->clone())); }
shared_ptr<DF_Full> DF_Full::copy() const { return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->copy())); }


void DF_Full::daxpy(const double a, const DF_Full& o) { data_->daxpy(a, o.data_); }
void DF_Full::daxpy(const double a, const DF_Half& o) { data_->daxpy(a, o.data_); }
void DF_Full::scale(const double a) { data_->scale(a); }
void DF_Full::symmetrize() { data_->symmetrize(); }


shared_ptr<DF_Full> DF_Full::apply_closed_2RDM() const {
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->apply_rhf_2RDM()));
}

// Caution
//   o strictly assuming that we are using natural orbitals.
shared_ptr<DF_Full> DF_Full::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, data_->apply_uhf_2RDM(amat, bmat)));
}

// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su
shared_ptr<DF_Half> DF_Full::back_transform(const double* c) const{
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc1_, data_->transform_third(c, df_->nbasis0(), true)));
}


shared_ptr<DF_AO> DF_Half::back_transform(const double* c) const{
  return shared_ptr<DF_AO>(new DF_AO(nbasis_, df_->nbasis1(), naux_, data_->transform_second(c, df_->nbasis1(), true)));
}


#endif
