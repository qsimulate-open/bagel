//
// Newint - Parallel electron correlation program.
// Filename: df.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <memory>
#include <src/util/f77.h>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;


void DensityFit::common_init(const vector<shared_ptr<Atom> >& atoms0,  const vector<vector<int> >& offsets0,
                             const vector<shared_ptr<Atom> >& atoms1,  const vector<vector<int> >& offsets1,
                             const vector<shared_ptr<Atom> >& aux_atoms,  const vector<vector<int> >& aux_offsets, const double throverlap, 
                             const bool compute_inverse) {

  // this will be distributed in the future.
  data_ = unique_ptr<double[]>(new double[nbasis0_*nbasis1_*naux_]);
  fill(data_.get(), data_.get()+nbasis0_*nbasis1_*naux_, 0.0);

  // info for basis 0
  vector<shared_ptr<Shell> > basis0;
  vector<int> offset0;
  int cnt = 0;
  for (auto aiter = atoms0.begin(); aiter != atoms0.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    basis0.insert(basis0.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = offsets0[cnt]; 
    offset0.insert(offset0.end(), tmpoff.begin(), tmpoff.end());
  }
  const int size0 = basis0.size();

  // info for basis 1
  vector<shared_ptr<Shell> > basis1;
  vector<int> offset1;
  cnt = 0;
  for (auto aiter = atoms1.begin(); aiter != atoms1.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    basis1.insert(basis1.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = offsets1[cnt]; 
    offset1.insert(offset1.end(), tmpoff.begin(), tmpoff.end());
  }
  const int size1 = basis1.size();

  // some info for auxiliary (i.e., DF) basis set
  vector<shared_ptr<Shell> > aux_basis; 
  vector<int> aux_offset;
  cnt = 0;
  for (auto aiter = aux_atoms.begin(); aiter != aux_atoms.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    aux_basis.insert(aux_basis.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = aux_offsets[cnt]; 
    aux_offset.insert(aux_offset.end(), tmpoff.begin(), tmpoff.end());
  }
  const int aux_size = aux_basis.size();

  if (basis0.front()->spherical() ^ basis1.front()->spherical()) throw runtime_error("do not mix spherical to cartesian...");
  const shared_ptr<Shell> b3(new Shell(basis0.front()->spherical()));

  if (basis0 == basis1) {
    assert(nbasis0_ == nbasis1_);
    for (int i0 = 0; i0 != size0; ++i0) {
      const shared_ptr<Shell>  b0 = basis0[i0];
      const int b0offset = offset0[i0]; 
      const int b0size = b0->nbasis();
      for (int i1 = i0; i1 != size1; ++i1) {
        const shared_ptr<Shell>  b1 = basis1[i1];
        const int b1offset = offset1[i1]; 
        const int b1size = b1->nbasis();
        for (int i2 = 0; i2 != aux_size; ++i2) {
          const shared_ptr<Shell>  b2 = aux_basis[i2];
          const int b2offset = aux_offset[i2]; 
          const int b2size = b2->nbasis();
  
          vector<shared_ptr<Shell> > input;
          input.push_back(b3);
          input.push_back(b2);
          input.push_back(b1);
          input.push_back(b0);
  
          // pointer to stack
          const double* ppt = compute_batch(input);
  
          // all slot in
          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2, ++ppt) {  
                data_[j2+naux_*(j1+nbasis1_*j0)] = data_[j2+naux_*(j0+nbasis1_*j1)] = *ppt;
              }
            }
          }
        }
      }
    }
  } else {
    assert(nbasis0_ == nbasis1_);
    for (int i0 = 0; i0 != size0; ++i0) {
      const shared_ptr<Shell>  b0 = basis0[i0];
      const int b0offset = offset0[i0]; 
      const int b0size = b0->nbasis();
      for (int i1 = 0; i1 != size1; ++i1) {
        const shared_ptr<Shell>  b1 = basis1[i1];
        const int b1offset = offset1[i1]; 
        const int b1size = b1->nbasis();
        for (int i2 = 0; i2 != aux_size; ++i2) {
          const shared_ptr<Shell>  b2 = aux_basis[i2];
          const int b2offset = aux_offset[i2]; 
          const int b2size = b2->nbasis();
  
          vector<shared_ptr<Shell> > input;
          input.push_back(b3);
          input.push_back(b2);
          input.push_back(b1);
          input.push_back(b0);
  
          // pointer to stack
          const double* ppt = compute_batch(input);
  
          // all slot in
          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2, ++ppt) {  
                data_[j2+naux_*(j1+nbasis1_*j0)] = *ppt;
              }
            }
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  data2_ = unique_ptr<double[]>(new double[naux_*naux_]);
  fill(data2(), data2()+naux_*naux_, 0.0);

  for (int i0 = 0; i0 != aux_size; ++i0) {
    const shared_ptr<Shell>  b0 = aux_basis[i0];
    const int b0offset = aux_offset[i0]; 
    const int b0size = b0->nbasis();
    for (int i1 = i0; i1 != aux_size; ++i1) {
      const shared_ptr<Shell>  b1 = aux_basis[i1];
      const int b1offset = aux_offset[i1]; 
      const int b1size = b1->nbasis();

      vector<shared_ptr<Shell> > input;
      input.push_back(b1);
      input.push_back(b3);
      input.push_back(b0);
      input.push_back(b3);

      // pointer to stack
      const double* ppt = compute_batch(input);
      for (int j0 = b0offset; j0 != b0offset+b0size; ++j0)
        for (int j1 = b1offset; j1 != b1offset+b1size; ++j1, ++ppt)
          data2_[j1+j0*naux_] = data2_[j0+j1*naux_] = *ppt;

    }
  }

  if (compute_inverse) {
    const size_t lwork = 5*naux_;
    unique_ptr<double[]> vec(new double[naux_]);
    unique_ptr<double[]> work(new double[max(lwork,naux_*naux_)]);

    int info;
    dsyev_("V", "U", naux_, data2_, naux_, vec, work, lwork, info); 
    if (info) throw runtime_error("dsyev failed in DF fock builder");

    for (int i = 0; i != naux_; ++i)
      vec[i] = vec[i] > throverlap ? 1.0/sqrt(sqrt(vec[i])) : 0.0;
    for (int i = 0; i != naux_; ++i) {
      for (int j = 0; j != naux_; ++j) {
        data2_[j+i*naux_] *= vec[i];
      }
    }
    // work now contains -1/2
    dgemm_("N", "T", naux_, naux_, naux_, 1.0, data2_, naux_, data2_, naux_, 0.0, work, naux_); 
    copy(work.get(), work.get()+naux_*naux_, data2());
  }

}

shared_ptr<DF_Half> DF_Half::clone() const {
  unique_ptr<double[]> dat(new double[size()]);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, dat));
}


shared_ptr<DF_Half> DF_Half::copy() const {
  unique_ptr<double[]> dat(new double[size()]);
  std::copy(data_.get(), data_.get()+size(), dat.get());
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, dat));
}


shared_ptr<DF_Half> DF_Half::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
  dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, d->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, out));
}


shared_ptr<DF_Half> DensityFit::compute_half_transform(const double* c, const size_t nocc) const {
  // it starts with DGEMM of inner index (for some reasons)...
  unique_ptr<double[]> tmp(new double[naux_*nbasis0_*nocc]);
  for (size_t i = 0; i != nbasis0_; ++i) {
    dgemm_("N", "N", naux_, nocc, nbasis1_, 1.0, data_.get()+i*naux_*nbasis1_, naux_, c, nbasis1_, 0.0, tmp.get()+i*naux_*nocc, naux_);
  }
  return shared_ptr<DF_Half>(new DF_Half(shared_from_this(), nocc, tmp));
}


void DF_Half::form_2index(unique_ptr<double[]>& target, const double a, const double b) const {
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, nbasis_, common, a, data_.get(), common, data_.get(), common, b, target.get(), nbasis_); 
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DF_Full> o, const double a, const double b) const {
  if (nocc_ != o->nocc1()) throw logic_error("nocc_ and o->nocc1() should be the same: in DF_Half::form_2index");
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, o->nocc2(), common, a, data_.get(), common, o->data(), common, b, target.get(), nbasis_); 
}


void DF_Half::form_4index(unique_ptr<double[]>& target) const {
  const int ndim = nbasis_ * nocc_;
  dgemm_("T", "N", ndim, ndim, naux_, 1.0, data_.get(), naux_, data_.get(), naux_, 0.0, target.get(), ndim); 
}

unique_ptr<double[]> DF_Half::form_4index() const {
  const size_t ndim = nbasis_ * nocc_;
  unique_ptr<double[]> out(new double[ndim*ndim]);
  form_4index(out);
  return move(out);
}

shared_ptr<DF_Full> DF_Half::compute_second_transform(const double* c, const size_t nocc) const {
  unique_ptr<double[]> tmp(new double[naux_*nocc_*nocc]);
  dgemm_("N", "N", naux_*nocc_, nocc, nbasis_, 1.0, data_.get(), naux_*nocc_, c, nbasis_, 0.0, tmp.get(), naux_*nocc_); 
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc_, nocc, tmp)); 
}


shared_ptr<DF_Full> DF_Full::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, d->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}

shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm) const {
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "T", naux_, nocc1_*nocc2_, nocc1_*nocc2_, 1.0, data_.get(), naux_, rdm, nocc1_*nocc2_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}



// forms all-internal 4-index MO integrals
void DF_Full::form_4index(unique_ptr<double[]>& target) const {
  const int dim = nocc1_ * nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, dim, naux, 1.0, data_.get(), naux, data_.get(), naux, 0.0, target.get(), dim); 
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
  dgemm_("T", "N", dim, odim, naux, 1.0, data_.get(), naux, o->data_.get(), naux, 0.0, target.get(), dim); 
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
  dgemm_("T", "N", dim, odim, naux, 1.0, data_.get(), naux, o->data_.get()+naux*odim*n, naux, 0.0, target.get(), dim); 
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
  dgemm_("T", "N", odim, dim, naux_, 1.0, o->data_3index(), naux_, tmp->data(), naux_, 0.0, target.get(), odim); 
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
  dgemm_("N", "T", naux_, naux_, nocc1_*nocc2_, 1.0, data(), naux_, o->data(), naux_, 0.0, out.get(), naux_); 
  return out;
}


shared_ptr<DF_Full> DF_Full::copy() const {
  unique_ptr<double[]> d(new double[size()]);
  std::copy(data(), data()+size(), d.get()); 
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


shared_ptr<DF_Full> DF_Full::apply_closed_2RDM() const {
  assert(nocc1_ == nocc2_);
  const int nocc = nocc1_; 
  unique_ptr<double[]> d(new double[size()]);
  fill(d.get(), d.get()+size(), 0.0);
  // exchange contributions
  daxpy_(size(), -2.0, data(), 1, d.get(), 1);
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, 1.0, data()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
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
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data(), naux_*nocc, amat, nocc, 0.0, d2.get(), naux_*nocc);
    for (int i = 0; i != nocc; ++i) {
      dgemm_("N", "N", naux_, nocc, nocc, -1.0, d2.get()+naux_*nocc*i, naux_, amat, nocc, 0.0, d.get()+naux_*nocc*i, naux_);
    }
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data(), naux_*nocc, bmat, nocc, 0.0, d2.get(), naux_*nocc);
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
    daxpy_(naux_, sum[i], data()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
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
  dgemm_("N", "T", naux_*nocc1_, nbas, nocc2_, 1.0, data(), naux_*nocc1_, c, nbas, 0.0, d.get(), naux_*nocc1_); 
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc1_, d)); 
}


shared_ptr<DF_AO> DF_Half::back_transform(const double* c) const{
  const int nbas = df_->nbasis0();
  unique_ptr<double[]> d(new double[nbas*nbasis_*naux_]);
  for (int i = 0; i != nbasis_; ++i)
    dgemm_("N", "T", naux_, nbas, nocc_, 1.0, data()+i*naux_*nocc_, naux_, c, nbas, 0.0, d.get()+i*naux_*nbas, naux_); 
  return shared_ptr<DF_AO>(new DF_AO(nbas, nbasis_, naux_, d)); 
}

