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


#include <src/util/f77.h>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>

using namespace std;

DensityFit::DensityFit(const int nbas, const int naux,
                                       const vector<shared_ptr<Atom> >& atoms,  const vector<vector<int> >& offsets,
                                       const vector<shared_ptr<Atom> >& aux_atoms,  const vector<vector<int> >& aux_offsets, const double throverlap)
   : nbasis_(nbas), naux_(naux) {

  // this will be distributed in the future.
  std::unique_ptr<double[]> data__(new double[nbasis_*nbasis_*naux]);
  data_ = move(data__);
  fill(data_.get(), data_.get()+nbasis_*nbasis_*naux, 0.0);

  vector<shared_ptr<Shell> > basis; 
  vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = offsets[cnt]; 
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }
  const int size = basis.size();

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

  const std::shared_ptr<Shell> b3(new Shell(basis.front()->spherical()));

  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<Shell>  b0 = basis[i0];
    const int b0offset = offset[i0]; 
    const int b0size = b0->nbasis();
    for (int i1 = i0; i1 != size; ++i1) {
      const shared_ptr<Shell>  b1 = basis[i1];
      const int b1offset = offset[i1]; 
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

        // TODO if I turn on primitive screening, it is broken.
        ERIBatch eribatch(input, 0.0);
        eribatch.compute();
        const double* eridata = eribatch.data();

        // all slot in
        for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
          for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
            for (int j2 = b2offset; j2 != b2offset + b2size; ++j2, ++eridata) {  
              data_[j2+naux*(j1+nbasis_*j0)] = data_[j2+naux*(j0+nbasis_*j1)] = *eridata;
            }
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  unique_ptr<double[]> data2__(new double[naux*naux]);
  data2_ = move(data2__);
  fill(data2(), data2()+naux*naux, 0.0);

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

      // TODO if I turn on primitive screening, it is broken.
      ERIBatch eribatch(input, 0.0);
      eribatch.compute();
      const double* eridata = eribatch.data();

      for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
        for (int j1 = b1offset; j1 != b1offset + b1size; ++j1, ++eridata) {  
          data2_[j1+j0*naux] = data2_[j0+j1*naux] = *eridata; 
        }
      }
    }
  }

  const int lwork = 5*naux;
  std::unique_ptr<double[]> vec(new double[naux]);
  std::unique_ptr<double[]> work(new double[std::max(lwork,naux*naux)]);

  int info;
  dsyev_("V", "U", naux, data2_, naux, vec, work, lwork, info); 
  if (info) throw std::runtime_error("dsyev failed in DF fock builder");

  for (int i = 0; i != naux; ++i)
    vec[i] = vec[i] > throverlap ? 1.0/std::sqrt(std::sqrt(vec[i])) : 0.0;
  for (int i = 0; i != naux; ++i) {
    for (int j = 0; j != naux; ++j) {
      data2_[j+i*naux] *= vec[i];
    }
  }
  // work now contains -1/2
  dgemm_("N", "T", naux, naux, naux, 1.0, data2_, naux, data2_, naux, 0.0, work, naux); 
  copy(work.get(), work.get()+naux*naux, data2());

}



shared_ptr<DF_Half> DensityFit::compute_half_transform(const double* c, const size_t nocc) const {
  unique_ptr<double[]> tmp(new double[naux_*nbasis_*nocc]);
  for (size_t i = 0; i != nbasis_; ++i) {
    dgemm_("N", "N", naux_, nocc, nbasis_, 1.0, data_.get()+i*naux_*nbasis_, naux_, c, nbasis_, 0.0, tmp.get()+i*naux_*nocc, naux_);
  }
  shared_ptr<DF_Half> out(new DF_Half(shared_from_this(), nocc, tmp));
  return out;
}


void DF_Half::form_2index(unique_ptr<double[]>& target, const double a, const double b) const {
  const int nbasis = df_->nbasis();
  const int common = nocc_ * df_->naux();
  dgemm_("T", "N", nbasis, nbasis, common, a, data_.get(), common, data_.get(), common, b, target.get(), nbasis); 
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DF_Full> o, const double a, const double b) const {
  assert(nocc_ == o->nocc1());
  const int nbasis = df_->nbasis();
  const int common = nocc_ * df_->naux();
  dgemm_("T", "N", nbasis, o->nocc2(), common, a, data_.get(), common, o->data(), common, b, target.get(), nbasis); 
}


void DF_Half::form_4index(unique_ptr<double[]>& target) const {
  const int ndim = df_->nbasis() * nocc_;
  const int naux = df_->naux();
  dgemm_("T", "N", ndim, ndim, naux, 1.0, data_.get(), naux, data_.get(), naux, 0.0, target.get(), ndim); 
}

unique_ptr<double[]> DF_Half::form_4index() const {
  const size_t ndim = df_->nbasis() * nocc_;
  unique_ptr<double[]> out(new double[ndim*ndim]);
  form_4index(out);
  return move(out);
}

shared_ptr<DF_Full> DF_Half::compute_second_transform(const double* c, const size_t nocc) const {
  const int naux = df_->naux();
  const int nbasis = df_->nbasis();
  unique_ptr<double[]> tmp(new double[naux*nocc_*nocc]);
  dgemm_("N", "N", naux*nocc_, nocc, nbasis, 1.0, data_.get(), naux*nocc_, c, nbasis, 0.0, tmp.get(), naux*nocc_); 
  shared_ptr<DF_Full> out(new DF_Full(df_, nocc_, nocc, tmp)); 
  return out;
}


void DF_Full::form_4index(unique_ptr<double[]>& target) const {
  const int dim = nocc1_ * nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, dim, naux, 1.0, data_.get(), naux, data_.get(), naux, 0.0, target.get(), dim); 
}


void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nocc1_ * o->nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, odim, naux, 1.0, data_.get(), naux, o->data_.get(), naux, 0.0, target.get(), dim); 
}


// Joperator. Note that (r,s) runs first; i.e., in the operator form
void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DensityFit> o) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nbasis() * o->nbasis();
  const int naux = df_->naux();
  shared_ptr<DF_Full> tmp = this->apply_J();
  dgemm_("T", "N", odim, dim, naux, 1.0, o->data_3index(), naux, tmp->data(), naux, 0.0, target.get(), odim); 
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DensityFit> o) const {
  const size_t dim = nocc1_ * o->nbasis();
  assert(nocc1_ == nocc2_);
  unique_ptr<double[]> out(new double[dim*dim]);
  form_4index(out, o);
  return move(out);
}

