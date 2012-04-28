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

using namespace std;


void DensityFit::common_init(const vector<shared_ptr<Atom> >& atoms0,  const vector<vector<int> >& offsets0,
                             const vector<shared_ptr<Atom> >& atoms1,  const vector<vector<int> >& offsets1,
                             const vector<shared_ptr<Atom> >& aux_atoms,  const vector<vector<int> >& aux_offsets, const double throverlap, const bool j2) {

  // this will be distributed in the future.
  std::unique_ptr<double[]> data__(new double[nbasis0_*nbasis1_*naux_]);
  data_ = move(data__);
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
  const std::shared_ptr<Shell> b3(new Shell(basis0.front()->spherical()));

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
  if (j2) {
    unique_ptr<double[]> data2__(new double[naux_*naux_]);
    data2_ = move(data2__);
    fill(data2(), data2()+naux_*naux_, 0.0);

    for (int i0 = 0; i0 != aux_size; ++i0) {
      const std::shared_ptr<Shell>  b0 = aux_basis[i0];
      const int b0offset = aux_offset[i0]; 
      const int b0size = b0->nbasis();
      for (int i1 = i0; i1 != aux_size; ++i1) {
        const std::shared_ptr<Shell>  b1 = aux_basis[i1];
        const int b1offset = aux_offset[i1]; 
        const int b1size = b1->nbasis();

        std::vector<std::shared_ptr<Shell> > input;
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

    const size_t lwork = 5*naux_;
    std::unique_ptr<double[]> vec(new double[naux_]);
    std::unique_ptr<double[]> work(new double[std::max(lwork,naux_*naux_)]);

    int info;
    dsyev_("V", "U", naux_, data2_, naux_, vec, work, lwork, info); 
    if (info) throw std::runtime_error("dsyev failed in DF fock builder");

    for (int i = 0; i != naux_; ++i)
      vec[i] = vec[i] > throverlap ? 1.0/std::sqrt(std::sqrt(vec[i])) : 0.0;
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



shared_ptr<DF_Half> DensityFit::compute_half_transform(const double* c, const size_t nocc) const {
  // it starts with DGEMM of inner index (for some reasons)...
  unique_ptr<double[]> tmp(new double[naux_*nbasis0_*nocc]);
  for (size_t i = 0; i != nbasis0_; ++i) {
    dgemm_("N", "N", naux_, nocc, nbasis1_, 1.0, data_.get()+i*naux_*nbasis1_, naux_, c, nbasis1_, 0.0, tmp.get()+i*naux_*nocc, naux_);
  }
  shared_ptr<DF_Half> out(new DF_Half(shared_from_this(), nocc, tmp));
  return out;
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
  shared_ptr<DF_Full> out(new DF_Full(df_, nocc_, nocc, tmp)); 
  return out;
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

