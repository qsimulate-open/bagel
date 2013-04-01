//
// BAGEL - Parallel electron correlation program.
// Filename: civec.cc
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


#include <src/fci/civec.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

Civec::Civec(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ = unique_ptr<double[]>(new double[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  fill_n(cc(), lena_*lenb_, 0.0);
}


Civec::Civec(shared_ptr<const Determinants> det, double* din_) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ptr_ = din_;
  fill_n(cc(), lena_*lenb_, 0.0);
}


Civec::Civec(const Civec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_) {
  cc_ = unique_ptr<double[]>(new double[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  copy_n(o.cc(), lena_*lenb_, cc());
}


// TODO Not efficient.
Civec::Civec(const DistCivec& o) : det_(o.det()), lena_(o.lena()), lenb_(o.lenb()) {
  cc_ = unique_ptr<double[]>(new double[size()]);
  cc_ptr_ = cc_.get();
  fill_n(cc_ptr_, size(), 0.0);
  copy_n(o.local(), o.asize()*lenb_, cc()+o.astart()*lenb_);
  mpi__->allreduce(cc_ptr_, size());
}


Civec::Civec(shared_ptr<Civec> o, shared_ptr<const Determinants> det) : det_(det), lena_(o->lena_), lenb_(o->lenb_) {
  assert(lena_ == det->lena() && lenb_ == det->lenb());
  cc_ = move(o->cc_);
  cc_ptr_ = cc_.get();
}


shared_ptr<Civec> Civec::transpose(shared_ptr<const Determinants> det) const {
  if (det == nullptr) det = det_->transpose();
  shared_ptr<Civec> ct(new Civec(det));
  double* cct = ct->data();
  mytranspose_(cc(), lenb_, lena_, cct);
  return ct;
}


double Civec::ddot(const Civec& other) const {
  assert( (lena_ == other.lena_) && (lenb_ == other.lenb_) );
  return ddot_(lena_*lenb_, cc(), 1, other.data(), 1);
}


void Civec::daxpy(double a, const Civec& other) {
  daxpy_(lena_*lenb_, a, other.data(), 1, cc(), 1);
}


double Civec::norm() const {
  return sqrt(ddot(*this));
}


void Civec::scale(const double a) {
  dscal_(lena_*lenb_, a, cc(), 1);
}


double Civec::variance() const {
  return ddot(*this) / (lena_*lenb_);
}


double Civec::orthog(list<shared_ptr<const Civec>> c) {
  for (auto& iter : c)
    project_out(iter);
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  scale(scal);
  return 1.0/scal;
}

double Civec::orthog(shared_ptr<const Civec> o) {
  list<shared_ptr<const Civec>> v = {o};
  return orthog(v);
}


Civec& Civec::operator/=(const Civec& o) {
  for (int i = 0; i != size(); ++i) {
    data(i) /= o.data(i);
  }
  return *this;
}

Civec Civec::operator/(const Civec& o) const {
  Civec out(*this);
  out /= o;
  return out;
}


shared_ptr<DistCivec> Civec::distcivec() const {
  shared_ptr<DistCivec> dist(new DistCivec(det_));
  copy_n(cc_ptr_+dist->astart()*lenb_, dist->asize()*lenb_, dist->local());
  return dist;
}

double Civec::spin_expectation() const {
  shared_ptr<Civec> S2 = spin();
  double out = ddot(*S2);

  return out;
}

// S^2 = S_z^2 + S_z + S_-S_+
shared_ptr<Civec> Civec::spin() const {
  shared_ptr<Civec> out(new Civec(det_));

  // First the easy part, S_z^2 + S_z
  const double sz = 0.5*static_cast<double>(det_->nspin());
  *out = *this;
  *out *= sz*sz + sz;

  const int norb = det_->norb();
  const int lena = det_->lena();
  const int lenb = det_->lenb();

  double* source = cc_ptr_;
  // This is a safe but probably slow implementation
  for (int aiter = 0; aiter < lena; ++aiter) {
    auto alphastring = det_->stringa(aiter);
    for (int biter = 0; biter < lenb; ++biter, ++source) {
      auto betastring = det_->stringb(biter);
      for (int i = 0; i < norb; ++i) {
        for (int j = 0; j < norb; ++j) {
          bitset<nbit__> abit = alphastring;
          bitset<nbit__> bbit = betastring;
          if (abit[i]) {
            abit.reset(i);
            if (!abit[j]) {
              abit.set(j);
              if (bbit[j]) {
                bbit.reset(j);
                if (!bbit[i]) {
                  bbit.set(i);

                  // Now the computation begins
                  const int atarget = det_->lexical<0>(abit);
                  const int btarget = det_->lexical<1>(bbit);
                  const int aphase = det_->sign(alphastring, j, i);
                  const int bphase = det_->sign(betastring, i, j);

                  out->element(btarget, atarget) -= static_cast<double>(aphase*bphase) * (*source);
                }
              }
            }
          }
        }

        if (betastring[i]) {
          out->element(biter,aiter) += *source;
        }
      }
    }
  }

  return out;
}

// S_- = \sum_i i_beta^\dagger i_alpha
shared_ptr<Civec> Civec::spin_lower(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = shared_ptr<Determinants>(new Determinants(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true));
  assert( (target_det->nelea() == det_->nelea()-1) && (target_det->neleb() == det_->neleb()+1) );
  shared_ptr<Civec> out(new Civec(target_det));

  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  double* source_data = cc_ptr_;
  // This is a safe but probably slow implementation
  for (int aiter = 0; aiter < source_lena; ++aiter) {
    auto alphastring = source_det->stringa(aiter);
    for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
      auto betastring = source_det->stringb(biter);
      for (int i = 0; i < norb; ++i) {
        bitset<nbit__> abit = alphastring;
        bitset<nbit__> bbit = betastring;
        if (abit[i]) {
          abit.reset(i);
          if (!bbit[i]) {
            bbit.set(i);

            const int atarget = target_det->lexical<0>(abit);
            const int btarget = target_det->lexical<1>(bbit);
                  // Now the computation begins

            const int aphase = source_det->sign(alphastring, i);
            const int bphase = source_det->sign(betastring, i);

            out->element(btarget, atarget) += static_cast<double>(aphase*bphase) * (*source_data);
          }
        }
      }
    }
  }

  return out;
}

// S_+ = \sum_i i_alpha^\dagger i_beta
shared_ptr<Civec> Civec::spin_raise(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = shared_ptr<Determinants>(new Determinants(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true));
  assert( (target_det->nelea() == det_->nelea()+1) && (target_det->neleb() == det_->neleb()-1) );
  shared_ptr<Civec> out(new Civec(target_det));

  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  double* source_data = cc_ptr_;
  // This is a safe but probably slow implementation
  for (int aiter = 0; aiter < source_lena; ++aiter) {
    auto alphastring = source_det->stringa(aiter);
    for (int biter = 0; biter < source_lenb; ++biter, ++source_data) {
      auto betastring = source_det->stringb(biter);
      for (int i = 0; i < norb; ++i) {
        bitset<nbit__> abit = alphastring;
        bitset<nbit__> bbit = betastring;
        if (bbit[i]) {
          bbit.reset(i);
          if (!abit[i]) {
            abit.set(i);

            const int atarget = target_det->lexical<0>(abit);
            const int btarget = target_det->lexical<1>(bbit);

            const int aphase = source_det->sign(alphastring, i);
            const int bphase = source_det->sign(betastring, i);

            out->element(btarget, atarget) += static_cast<double>(aphase*bphase) * (*source_data);
          }
        }
      }
    }
  }

  return out;
}
