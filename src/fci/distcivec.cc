//
// BAGEL - Parallel electron correlation program.
// Filename: distcivec.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/fci/civec.h>
#include <src/math/algo.h>
#include <src/parallel/distqueue.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

namespace bagel { namespace DFCI {
  class SpinTask {
    protected:
      const bitset<nbit__> abit_;
      double* target_;
      const DistCivector<double>* cc_;

      unique_ptr<double[]> buf_;
      vector<int> requests_;

    public:
      SpinTask(const bitset<nbit__> a, double* t, const DistCivector<double>* c) : abit_(a), target_(t), cc_(c) {
        shared_ptr<const Determinants> det = c->det();
        const size_t lb = det->lenb();
        const int norb = det->norb();

        buf_ = unique_ptr<double[]>(new double[lb * abit_.count() * (norb - abit_.count() + 1)]);

        for (int i = 0, k = 0; i < norb; ++i) {
          for (int j = 0; j < norb; ++j) {
            if ( abit_[j] ) {
              bitset<nbit__> tarbit = abit_; tarbit.reset(j);
              if ( !tarbit[i] ) {
                tarbit.set(i);
                const int l = cc_->get_bstring_buf(buf_.get() + lb*k++, det->lexical<0>(tarbit));
                if (l >= 0) requests_.push_back(l);
              }
            }
          }
        }
      }

      bool test() {
        bool out = true;
        for (auto i = requests_.begin(); i != requests_.end(); ) {
          if (mpi__->test(*i)) {
            i = requests_.erase(i);
          }
          else {
            ++i;
            out = false;
          }
        }
        return out;
      }

      void compute() {
        shared_ptr<const Determinants> det = cc_->det();

        const size_t lb = det->lenb();
        const int norb = det->norb();

        for (int i = 0, k = 0; i < norb; ++i) {
          for (int j = 0; j < norb; ++j) {
            if ( abit_[j] ) {
              bitset<nbit__> tarbit = abit_; tarbit.reset(j);
              if ( !tarbit[i] ) {
                const double* source = buf_.get() + lb*k++;
                const int aphase = det->sign(abit_, i, j);
                for (auto& iter : det->phib(j, i))
                  target_[iter.target] -= static_cast<double>(aphase * iter.sign) * source[iter.source];
              }
            }
          }
        }
      }
  };
} }

template <>
shared_ptr<DistCivector<double>> DistCivector<double>::spin() const {
  const size_t la = det_->lena();
  const size_t lb = det_->lenb();
  auto out = this->clone();

  this->init_mpi_recv();

#ifndef USE_SERVER_THREAD
  DistQueue<DFCI::SpinTask, const DistCivector<double>*> dq(this);
#else
  DistQueue<DFCI::SpinTask> dq;
#endif

  for (size_t ia = astart_; ia < aend_; ++ia)
    dq.emplace_and_compute(det_->stringa(ia), out->local() + (ia - astart_) * lb, this);

  const double sz = 0.5*static_cast<double>(det_->nspin());
  out->ax_plus_y(sz*sz + sz + static_cast<double>(det_->neleb()), *this);

  dq.finish();

  this->terminate_mpi_recv();
  return out;
}

template<>
void DistCivector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<DistCivector<double>> S2 = spin();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while( fabs(actual_expectation - pure_expectation) > thresh ) {
    if ( k > max_spin ) { this->print(0.05); throw runtime_error("Spin decontamination failed."); }

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2); 

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    actual_expectation = dot_product(*S2);

    k += 2;
  }
}

namespace bagel { namespace DFCI {
  class LowerTask {
    protected:
      const bitset<nbit__> abit_;
      double* target_;
      const DistCivector<double>* cc_;
      shared_ptr<const Determinants> tdet_;

      unique_ptr<double[]> buf_;
      vector<int> requests_;

    public:
      LowerTask(const bitset<nbit__> a, double* t, const DistCivector<double>* c, shared_ptr<const Determinants> td) : abit_(a), target_(t), cc_(c), tdet_(td) {
        shared_ptr<const Determinants> det = c->det();
        const size_t lbs = det->lenb();
        const int norb = det->norb();

        buf_ = unique_ptr<double[]>(new double[lbs * (norb - abit_.count())]);

        for (int i = 0, k = 0; i < norb; ++i) {
          if ( !abit_[i] ) {
            bitset<nbit__> sourcebit = abit_; sourcebit.set(i);
            const int l = cc_->get_bstring_buf(buf_.get() + lbs*k++, det->lexical<0>(sourcebit));
            if (l >= 0) requests_.push_back(l);
          }
        }
      }

      bool test() {
        bool out = true;
        for (auto i = requests_.begin(); i != requests_.end(); ) {
          if (mpi__->test(*i)) {
            i = requests_.erase(i);
          }
          else {
            ++i;
            out = false;
          }
        }
        return out;
      }

      void compute() {
        shared_ptr<const Determinants> sdet = cc_->det();

        const size_t lbs = sdet->lenb();
        const int norb = sdet->norb();

        for (int i = 0, k = 0; i < norb; ++i) {
          if ( abit_[i] ) continue;
          bitset<nbit__> sabit = abit_; sabit.set(i);
          const int aphase = sdet->sign<0>(sabit, i);
          const double* sourcedata = buf_.get() + lbs * k++;
          double* targetdata = target_;
          for ( auto& bbit : tdet_->stringb() ) {
            if ( bbit[i] ) {
              bitset<nbit__> sbbit = bbit; sbbit.reset(i);
              const int bphase = -sdet->sign<1>(sbbit, i);
              *targetdata += static_cast<double>(aphase * bphase) * sourcedata[sdet->lexical<1>(sbbit)];
            }
            ++targetdata;
          }
        }
      }
  };
} }

template<>
shared_ptr<DistCivector<double>> DistCivector<double>::spin_lower(shared_ptr<const Determinants> tdet) const {
  if (!tdet) tdet = make_shared<Determinants>(det_->norb(), det_->nelea()-1, det_->neleb()+1, det_->compress(), true);
  assert( (tdet->nelea() == det_->nelea()-1) && (tdet->neleb() == det_->neleb()+1) );

  auto out = make_shared<DistCivector<double>>(tdet);
  this->init_mpi_recv();

  const size_t tastart = out->astart();
  const size_t taend = out->aend();
  const size_t lbt = tdet->lenb();

#ifndef USE_SERVER_THREAD
  DistQueue<DFCI::LowerTask, const DistCivector<double>*> dq(this);
#else
  DistQueue<DFCI::LowerTask> dq;
#endif

  for (size_t ia = tastart; ia < taend; ++ia)
    dq.emplace_and_compute(tdet->stringa(ia), out->local() + (ia - tastart) * lbt, this, tdet);

  dq.finish();

  this->terminate_mpi_recv();
  return out;
}

namespace bagel { namespace DFCI {
  class RaiseTask {
    protected:
      const bitset<nbit__> abit_;
      double* target_;
      const DistCivector<double>* cc_;
      shared_ptr<const Determinants> tdet_;

      unique_ptr<double[]> buf_;
      vector<int> requests_;

    public:
      RaiseTask(const bitset<nbit__> a, double* t, const DistCivector<double>* c, shared_ptr<const Determinants> td) : abit_(a), target_(t), cc_(c), tdet_(td) {
        shared_ptr<const Determinants> det = c->det();
        const size_t lbs = det->lenb();
        const int norb = det->norb();

        buf_ = unique_ptr<double[]>(new double[lbs * abit_.count()]);

        for (int i = 0, k = 0; i < norb; ++i) {
          if ( abit_[i] ) {
            bitset<nbit__> sourcebit = abit_; sourcebit.reset(i);
            const int l = cc_->get_bstring_buf(buf_.get() + lbs*k++, det->lexical<0>(sourcebit));
            if (l >= 0) requests_.push_back(l);
          }
        }
      }

      bool test() {
        bool out = true;
        for (auto i = requests_.begin(); i != requests_.end(); ) {
          if (mpi__->test(*i)) {
            i = requests_.erase(i);
          }
          else {
            ++i;
            out = false;
          }
        }
        return out;
      }

      void compute() {
        shared_ptr<const Determinants> sdet = cc_->det();

        const size_t lbs = sdet->lenb();
        const int norb = sdet->norb();

        for (int i = 0, k = 0; i < norb; ++i) {
          if ( !abit_[i] ) continue;
          bitset<nbit__> sabit = abit_; sabit.reset(i);
          const int aphase = sdet->sign<0>(sabit, i);
          const double* sourcedata = buf_.get() + lbs * k++;
          double* targetdata = target_;
          for ( auto& bbit : tdet_->stringb() ) {
            if ( !bbit[i] ) {
              bitset<nbit__> sbbit = bbit; sbbit.set(i);
              const int bphase = sdet->sign<1>(sbbit, i);
              *targetdata += static_cast<double>(aphase * bphase) * sourcedata[sdet->lexical<1>(sbbit)];
            }
            ++targetdata;
          }
        }
      }
  };
} }

template<>
shared_ptr<DistCivector<double>> DistCivector<double>::spin_raise(shared_ptr<const Determinants> tdet) const {
  if (!tdet) tdet = make_shared<Determinants>(det_->norb(), det_->nelea()+1, det_->neleb()-1, det_->compress(), true);
  assert( (tdet->nelea() == det_->nelea()+1) && (tdet->neleb() == det_->neleb()-1) );

  auto out = make_shared<DistCivector<double>>(tdet);
  this->init_mpi_recv();

  const size_t tastart = out->astart();
  const size_t taend = out->aend();
  const size_t lbt = out->lenb();

#ifndef USE_SERVER_THREAD
  DistQueue<DFCI::RaiseTask, const DistCivector<double>*> dq(this);
#else
  DistQueue<DFCI::RaiseTask> dq;
#endif

  for (size_t ia = tastart; ia < taend; ++ia)
    dq.emplace_and_compute(tdet->stringa(ia), out->local() + (ia - tastart) * lbt, this, tdet);

  dq.finish();

  this->terminate_mpi_recv();
  return out;
}
