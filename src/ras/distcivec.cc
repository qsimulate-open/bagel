//
// BAGEL - Parallel electron correlation program.
// Filename: ras/distcivector.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <iomanip>
#include <unordered_map>

#include <src/ras/distcivector.h>
#include <src/parallel/distqueue.h>

using namespace std;
using namespace bagel;

// Computes <S^2>

namespace bagel {
  namespace RAS {
    struct DistSpinTask {
      const std::bitset<nbit__> abit_;
      const DistRASCivector<double>* this_;
      shared_ptr<DistRASCivector<double>> out_;
      shared_ptr<const RASDeterminants> det_;
      unordered_map<size_t, size_t>* lexicalmap_;

      unique_ptr<double[]> buf_;
      list<int> requests_;

      DistSpinTask(const std::bitset<nbit__> t, const DistRASCivector<double>* th, shared_ptr<DistRASCivector<double>> o, shared_ptr<const RASDeterminants> d, unordered_map<size_t, size_t>* lex) :
        abit_(t), this_(th), out_(o), det_(d), lexicalmap_(lex)
      {
        const size_t lb = det_->lenb();
        const int norb = det_->norb();

        const size_t alexical = det_->lexical<0>(abit_);
        buf_ = unique_ptr<double[]>(new double[lb * det_->phia(alexical).size()]);

        int k = 0;
        for (auto& iter : det_->phia(alexical)) {
          // loop through blocks
          for (auto& iblock : this_->allowed_blocks<0>(det_->stringa(iter.source))) {
            const int l = iblock->get_bstring_buf(buf_.get() + lb*k + iblock->stringb()->offset(), iter.source - iblock->stringa()->offset());
            if (l >= 0) requests_.push_back(l);
          }
          ++k;
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
        const int norb = det_->norb();
        const size_t lex_a = det_->lexical<0, 0>(abit_);

        vector<shared_ptr<DistRASBlock<double>>> allowed_blocks = out_->allowed_blocks<0>(abit_);

        int k = 0;
        for (auto& iter : det_->phia(det_->lexical<0>(abit_))) {
          const int j = iter.ij / norb;
          const int i = iter.ij % norb;
          bitset<nbit__> mask1; mask1.set(j); mask1.set(i);
          bitset<nbit__> mask2; mask2.set(j);

          bitset<nbit__> maskij; maskij.set(j); maskij.flip(i);

          const double* source = buf_.get() + det_->lenb() * k++;
          for (auto& iblock : allowed_blocks) {
            const size_t lb = iblock->lenb();
            double* odata = iblock->local() + lb * (lex_a - iblock->astart());
            for (auto& ib : *iblock->stringb()) {
              if ( ((ib & mask1) ^ mask2).none() ) { // equivalent to "ib[j] && (ii == jj || !ib[i])"
                const bitset<nbit__> bsostring = ib ^ maskij;
                *odata -= static_cast<double>(iter.sign * det_->sign(bsostring, i, j)) * source[(*lexicalmap_)[bsostring.to_ullong()]];
              }
              ++odata;
            }
          }
        }
      }
    };
  }
}

// Returns S^2 | civec >
// S^2 = S_z^2 + S_z + S_-S_+ with S_-S_+ = nbeta - \sum_{ij} j_alpha^dagger i_alpha i_beta^dagger j_beta
template<>
shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin() const {
  auto out = make_shared<DistRASCivector<double>>(det_);

  unordered_map<size_t, size_t> lexicalmap;
  for (size_t i = 0; i < det_->lenb(); ++i)
    lexicalmap[det_->stringb(i).to_ullong()] = i;

  this->init_mpi_recv();

#ifndef USE_SERVER_THREAD
    DistQueue<RAS::DistSpinTask, const DistRASCivector<double>*> tasks(this);
#else
    DistQueue<RAS::DistSpinTask> tasks;
#endif

  const size_t nthreads = resources__->max_num_threads();
  auto spin_task = [&nthreads, &out, &lexicalmap, &tasks, this] (const size_t myid) {

    for (auto& ispace : this->det()->stringspacea()) {
      if (!ispace) continue;
      StaticDist global_dist(ispace->size(), mpi__->size());
      size_t astart, aend;
      tie(astart, aend) = global_dist.range(mpi__->rank());

      StaticDist thread_dist(aend - astart, nthreads);
      size_t thread_start, thread_end;
      tie(thread_start, thread_end) = thread_dist.range(myid);

      for (size_t ia = thread_start; ia < thread_end; ++ia)
        tasks.emplace_and_compute(this->det()->stringa(ia + astart + ispace->offset()), this, out, this->det(), &lexicalmap);
    }

    ( myid == 0 ) ? tasks.finish_master() : tasks.finish_worker();
  };

  vector<thread> thread_tasks;
  for (size_t ith = 1; ith < nthreads; ++ith)
    thread_tasks.push_back( thread(spin_task, ith) );

  // This is master thread
  spin_task(0);

  const double sz = static_cast<double>(det_->nspin()) * 0.5;
  const double fac = sz*sz + sz + static_cast<double>(det_->neleb());

  out->ax_plus_y(fac, *this);

  for (auto& th : thread_tasks)
    th.join();

  this->terminate_mpi_recv();

  return out;
}

// S_- = \sum_i i^dagger_beta i_alpha
template<> shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
#if 0
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert( (tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1) );
  auto out = make_shared<DistRASCivec>(tdet);

  const int norb = sdet->norb();
  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // maps bits to their local offsets
  unordered_map<size_t, size_t> alex;
  for (auto& ispace : sdet->stringspacea()) {
    if (ispace)
      for (auto& abit : *ispace) alex[abit.to_ullong()] = ispace->lexical<0>(abit);
  }

  unordered_map<size_t, size_t> blex;
  for (auto& ispace : sdet->stringspaceb()) {
    if (ispace)
      for (auto& bbit : *ispace) blex[bbit.to_ullong()] = ispace->lexical<0>(bbit);
  }

  auto lower_ras = [&sdet, &alex, &blex] (shared_ptr<const RASBlock<double>> sblock, shared_ptr<RASBlock<double>> tblock, const int nstart, const int nfence) {
    const size_t lb = sblock->lenb();
    double* odata = tblock->data();
    for (auto& abit : *tblock->stringa()) {
      for (auto& bbit : *tblock->stringb()) {
        for ( int i = nstart; i < nfence; ++i) {
          if (abit[i] || !bbit[i]) continue;
          bitset<nbit__> sabit = abit; sabit.set(i);
          bitset<nbit__> sbbit = bbit; sbbit.reset(i);

          const double phase = static_cast<double>(-1 * sdet->sign<0>(sabit, i) * sdet->sign<1>(sbbit,i));

          *odata += phase * sblock->element( blex[sbbit.to_ullong()] + alex[sabit.to_ullong()] * lb );
        }
        ++odata;
      }
    }
  };

  // The important thing to notice is that for all orbitals in a single RAS space, each block in the source is sent to a single block in target
  for (auto& iblock : out->blocks()) {
    if (!iblock) continue;
    const int nha = iblock->stringa()->nholes();
    const int nhb = iblock->stringb()->nholes();
    const int npa = iblock->stringa()->nparticles();
    const int npb = iblock->stringb()->nparticles();
    const int n2a = iblock->stringa()->nele2();
    const int n2b = iblock->stringb()->nele2();

    if ( (ras1 > 0) && (nhb < ras1) && (nha > 0) ) lower_ras(this->block(nha-1,nhb+1,npa,npb), iblock, 0, ras1);
    if ( (ras2 > 0) && (n2b > 0) && (n2a < ras2) ) lower_ras(this->block(nha, nhb, npa, npb), iblock, ras1, ras1 + ras2);
    if ( (ras3 > 0) && (npb > 0) && (npa < ras3) ) lower_ras(this->block(nha, nhb, npa+1, npb-1), iblock, ras1+ras2, ras1+ras2+ras3);
  }

  return out;
#else
  return shared_ptr<DistRASCivector<double>>();
#endif
}

// S_+ = \sum_i i^dagger_alpha i_beta
template<> shared_ptr<DistRASCivector<double>> DistRASCivector<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
#if 0
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert( (tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1) );
  auto out = make_shared<DistRASCivec>(tdet);

  const int norb = sdet->norb();
  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // maps bits to their local offsets
  unordered_map<size_t, size_t> alex;
  for (auto& ispace : det_->stringspacea()) {
    if (ispace)
      for (auto& abit : *ispace) alex[abit.to_ullong()] = ispace->lexical<0>(abit);
  }

  unordered_map<size_t, size_t> blex;
  for (auto& ispace : det_->stringspaceb()) {
    if (ispace)
      for (auto& bbit : *ispace) blex[bbit.to_ullong()] = ispace->lexical<0>(bbit);
  }

  auto raise_ras = [&sdet, &alex, &blex] (shared_ptr<const RASBlock<double>> sblock, shared_ptr<RASBlock<double>> tblock, const int nstart, const int nfence) {
    const size_t lb = sblock->lenb();
    double* odata = tblock->data();
    for (auto& abit : *tblock->stringa()) {
      for (auto& bbit : *tblock->stringb()) {
        for ( int i = nstart; i < nfence; ++i) {
          if (!abit[i] || bbit[i]) continue;
          bitset<nbit__> sabit = abit; sabit.reset(i);
          bitset<nbit__> sbbit = bbit; sbbit.set(i);

          const double phase = static_cast<double>(sdet->sign<0>(sabit, i) * sdet->sign<1>(sbbit,i));

          *odata += phase * sblock->element( blex[sbbit.to_ullong()] + alex[sabit.to_ullong()] * lb );
        }
        ++odata;
      }
    }
  };

  // The important thing to notice is that for all orbitals in a single RAS space, each block in the source is sent to a single block in target
  for (auto& iblock : out->blocks()) {
    if (!iblock) continue;
    const int nha = iblock->stringa()->nholes();
    const int nhb = iblock->stringb()->nholes();
    const int npa = iblock->stringa()->nparticles();
    const int npb = iblock->stringb()->nparticles();
    const int n2a = iblock->stringa()->nele2();
    const int n2b = iblock->stringb()->nele2();

    if ( (ras1 > 0) && (nha < ras1) && (nhb > 0) ) raise_ras(this->block(nha+1,nhb-1,npa,npb), iblock, 0, ras1);
    if ( (ras2 > 0) && (n2a > 0) && (n2b < ras2) ) raise_ras(this->block(nha, nhb, npa, npb), iblock, ras1, ras1 + ras2);
    if ( (ras3 > 0) && (npa > 0) && (npb < ras3) ) raise_ras(this->block(nha, nhb, npa-1, npb+1), iblock, ras1+ras2, ras1+ras2+ras3);
  }

  return out;
#else
  return shared_ptr<DistRASCivector<double>>();
#endif
}

template<> void DistRASCivector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<DistRASCivec> S2 = spin();
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
