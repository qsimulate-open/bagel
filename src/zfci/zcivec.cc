//
// BAGEL - Parallel electron correlation program.
// Filename: zcivec.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
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


#include <src/zfci/zcivec.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


ZCivec::ZCivec(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ = unique_ptr<complex<double>[]>(new complex<double>[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  fill_n(cc(), lena_*lenb_, 0.0);
}


ZCivec::ZCivec(shared_ptr<const Determinants> det, complex<double>* din_) : det_(det), lena_(det->lena()), lenb_(det->lenb()) {
  cc_ptr_ = din_;
  fill_n(cc(), lena_*lenb_, 0.0);
}


ZCivec::ZCivec(const ZCivec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_) {
  cc_ = unique_ptr<complex<double>[]>(new complex<double>[lena_*lenb_]);
  cc_ptr_ = cc_.get();
  copy_n(o.cc(), lena_*lenb_, cc());
}


#ifndef ZDISTCIVEC_NOT_IMPLELEMTED
// TODO Not efficient.
ZCivec::ZCivec(const ZDistCivec& o) : det_(o.det()), lena_(o.lena()), lenb_(o.lenb()) {
  cc_ = unique_ptr<complex<double>[]>(new complex<double>[size()]);
  cc_ptr_ = cc_.get();
  fill_n(cc_ptr_, size(), 0.0);
  copy_n(o.local(), o.asize()*lenb_, cc()+o.astart()*lenb_);
  mpi__->allreduce(cc_ptr_, size());
}
#endif


ZCivec::ZCivec(shared_ptr<ZCivec> o, shared_ptr<const Determinants> det) : det_(det), lena_(o->lena_), lenb_(o->lenb_) {
  assert(lena_ == det->lena() && lenb_ == det->lenb());
  cc_ = move(o->cc_);
  cc_ptr_ = cc_.get();
}


shared_ptr<ZCivec> ZCivec::transpose(shared_ptr<const Determinants> det) const {
  if (det == nullptr) det = det_->transpose();
  auto ct = make_shared<ZCivec>(det);
  complex<double>* cct = ct->data();
  mytranspose_complex_(cc(), lenb_, lena_, cct);
  return ct;
}


complex<double> ZCivec::zdotc(const ZCivec& other) const {
  assert( (lena_ == other.lena_) && (lenb_ == other.lenb_) );
  return zdotc_(lena_*lenb_, cc(), 1, other.data(), 1);
}


void ZCivec::zaxpy(complex<double> a, const ZCivec& other) {
  zaxpy_(lena_*lenb_, a, other.data(), 1, cc(), 1);
}


double ZCivec::norm() const {
  complex<double> normc = zdotc(*this);
  assert(abs(normc.imag()) < 1.0e-8);
  return sqrt(normc.real());
}


void ZCivec::scale(const complex<double> a) {
  zscal_(lena_*lenb_, a, cc(), 1);
}


double ZCivec::variance() const {
  const double n = norm();
  return n*n / (lena_*lenb_);
}


complex<double> ZCivec::orthog(list<shared_ptr<const ZCivec>> c) {
  for (auto& iter : c)
    project_out(iter);
  const double norm = this->norm();
  const double scal = (abs(norm*norm)<1.0e-60 ? 0.0 : 1.0/norm);
  scale(scal);
  return 1.0/scal;
}

complex<double> ZCivec::orthog(shared_ptr<const ZCivec> o) {
  list<shared_ptr<const ZCivec>> v = {o};
  return orthog(v);
}


ZCivec& ZCivec::operator/=(const ZCivec& o) {
  for (int i = 0; i != size(); ++i) {
    data(i) /= o.data(i);
  }
  return *this;
}

ZCivec ZCivec::operator/(const ZCivec& o) const {
  ZCivec out(*this);
  out /= o;
  return out;
}


#ifndef ZDISTCIVEC_NOT_IMPLELEMTED
shared_ptr<ZDistCivec> ZCivec::zdistcivec() const {
  auto dist = make_shared<ZDistCivec>(det_);
  copy_n(cc_ptr_+dist->astart()*lenb_, dist->asize()*lenb_, dist->local());
  return dist;
}
#endif

double ZCivec::spin_expectation() const {
  shared_ptr<ZCivec> S2 = spin();
  complex<double> out = zdotc(*S2);
  assert(abs(out.imag()) < 1e-8);
  return out.real();
}


void ZCivec::print(const double thr) const {
  const complex<double>* i = cc();
  // multimap sorts elements so that they will be shown in the descending order in magnitude
  multimap<double, tuple<complex<double>, bitset<nbit__>, bitset<nbit__>>> tmp;
  for (auto& ia : det_->stringa()) {
    for (auto& ib : det_->stringb()) {
      if (abs(*i) > thr) {
        tmp.insert(make_pair(-sqrt(pow((*i).real(),2)+pow((*i).imag(),2)), make_tuple(*i, ia, ib)));
      }
      ++i;
    }
  }
  for (auto& iter : tmp) {
    cout << "       " << det_->print_bit(get<1>(iter.second), get<2>(iter.second))
         << "  " << setprecision(10) << setw(15) << get<0>(iter.second) << endl;
  }
}

//TODO check that this is correct
// S^2 = S_z^2 + S_z + S_-S_+
shared_ptr<ZCivec> ZCivec::spin() const {
  auto out = make_shared<ZCivec>(det_);

  // First the easy part, S_z^2 + S_z
  const double sz = 0.5*static_cast<double>(det_->nspin());
  *out = *this;
  *out *= sz*sz + sz + det_->neleb();

  const int norb = det_->norb();
  const int lena = det_->lena();
  const int lenb = det_->lenb();

  auto intermediate = make_shared<ZCivec>(det_);

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      intermediate->zero();

      for ( auto& iter : det_->phia(i,j) ) {
        const complex<double>* source = this->element_ptr(0, iter.source);
        complex<double>* target = intermediate->element_ptr(0, iter.target);
        double sign = static_cast<double>(iter.sign);
        zaxpy_(lenb, sign, source, 1, target, 1);
      }
      for ( auto& iter : det_->phib(j,i) ) {
        complex<double>* source = intermediate->element_ptr(iter.source, 0);
        complex<double>* target = out->element_ptr(iter.target, 0);
        double sign = static_cast<double>(-iter.sign);
        zaxpy_(lena, sign, source, lenb, target, lenb);
      }
    }
  }
  return out;
}

// S_- = \sum_i i_beta^\dagger i_alpha
shared_ptr<ZCivec> ZCivec::spin_lower(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);
  assert( (target_det->nelea() == det_->nelea()-1) && (target_det->neleb() == det_->neleb()+1) );
  auto out = make_shared<ZCivec>(target_det);

  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  complex<double>* source_data = cc_ptr_;
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

            const int aphase = source_det->sign<0>(alphastring, i);
            const int bphase = -1*source_det->sign<1>(betastring, i);

            out->element(btarget, atarget) += static_cast<complex<double>>(aphase*bphase) * (*source_data);
          }
        }
      }
    }
  }

  return out;
}

// S_+ = \sum_i i_alpha^\dagger i_beta
shared_ptr<ZCivec> ZCivec::spin_raise(shared_ptr<const Determinants> target_det) const {
  if (target_det == nullptr)
    target_det = make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);
  assert( (target_det->nelea() == det_->nelea()+1) && (target_det->neleb() == det_->neleb()-1) );
  auto out = make_shared<ZCivec>(target_det);

  shared_ptr<const Determinants> source_det = det_;
  const int norb = source_det->norb();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  complex<double>* source_data = cc_ptr_;
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

            const int aphase = source_det->sign<0>(alphastring, i);
            const int bphase = source_det->sign<1>(betastring, i);

            out->element(btarget, atarget) += static_cast<complex<double>>(aphase*bphase) * (*source_data);
          }
        }
      }
    }
  }

  return out;
}

//TODO check this
void ZCivec::spin_decontaminate(const double thresh) {

  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<ZCivec> S2 = spin();

  int k = nspin + 2;
  while( abs(zdotc(*S2) - expectation) > thresh ) {
    if ( k > max_spin ) throw runtime_error("Spin decontamination failed.");

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    zaxpy(factor, *S2);

    const double norm = this->norm();
    const double rescale = (abs(norm*norm) > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();

    k += 2;
  }
}
