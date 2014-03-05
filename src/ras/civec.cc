//
// BAGEL - Parallel electron correlation program.
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/ras/civector.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;

// Computes <S^2>
template<>
double RASCivector<double>::spin_expectation() const {
  shared_ptr<const RASCivector<double>> S2 = spin();
  return dot_product(*S2);
}

namespace bagel {
  namespace RAS {
    struct SpinTask {
      const std::bitset<nbit__> target_;
      const RASCivector<double>* this_;
      shared_ptr<RASCivector<double>> out_;
      shared_ptr<const RASDeterminants> det_;
      unordered_map<bitset<nbit__>, size_t>* lexicalmap_;

      SpinTask(const std::bitset<nbit__> t, const RASCivector<double>* th, shared_ptr<RASCivector<double>> o, shared_ptr<const RASDeterminants> d, unordered_map<bitset<nbit__>, size_t>* lex) :
        target_(t), this_(th), out_(o), det_(d), lexicalmap_(lex) {}

      void compute() {
        const int norb = det_->norb();

        unique_ptr<double[]> source(new double[det_->lenb()]);

        for (auto& iter : det_->phia(det_->lexical_offset<0>(target_))) {
          const int ii = iter.ij / norb;
          const int jj = iter.ij % norb;
          bitset<nbit__> mask1; mask1.set(ii); mask1.set(jj);
          bitset<nbit__> mask2; mask2.set(ii);

          bitset<nbit__> maskij; maskij.set(ii); maskij.flip(jj);

          fill_n(source.get(), det_->lenb(), 0.0);
          vector<shared_ptr<RASBlock<double>>> sourceblocks = out_->allowed_blocks<0>(det_->string_bits_a(iter.source));
          for (auto& iblock : sourceblocks) {
            const size_t offset = iblock->stringsb()->offset();
            copy_n(&this_->element(iblock->stringsb()->strings(0), det_->string_bits_a(iter.source)), iblock->lenb(), source.get()+offset);
          }

          for (auto& iblock : out_->allowed_blocks<0>(target_)) {
            double* outelement = &out_->element(iblock->stringsb()->strings(0), target_);
            for (auto& btstring : *iblock->stringsb()) {
              if ( ((btstring & mask1) ^ mask2).none() ) { // equivalent to "btstring[ii] && (ii == jj || !btstring[jj])"
                const bitset<nbit__> bsostring = btstring ^ maskij;
                if (det_->allowed(det_->string_bits_a(iter.source), bsostring))
                  *outelement -= static_cast<double>(iter.sign * det_->sign(bsostring, ii, jj)) * source[(*lexicalmap_)[bsostring]];
              }
              ++outelement;
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
shared_ptr<RASCivector<double>> RASCivector<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);

  unordered_map<bitset<nbit__>, size_t> lexicalmap;
  for (auto& i : det_->string_bits_b())
    lexicalmap[i] = det_->lexical_offset<1>(i);

  TaskQueue<RAS::SpinTask> tasks(det_->string_bits_a().size());

  for (auto& istring : det_->string_bits_a()) {
    tasks.emplace_back(istring, this, out, det_, &lexicalmap);
  }

  tasks.compute();

  const double sz = static_cast<double>(det_->nspin()) * 0.5;
  const double fac = sz*sz + sz + static_cast<double>(det_->neleb());

  out->ax_plus_y(fac, *this);

  return out;
}

// S_- = \sum_i i^dagger_beta i_alpha
template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert( (tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1) );
  auto out = make_shared<RASCivec>(tdet);

  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // maps bits to their local offsets
  unordered_map<bitset<nbit__>, size_t> alex;
  for (auto& ispace : *sdet->stringspacea()) {
    for (auto& abit : *ispace) alex[abit] = ispace->lexical_zero(abit);
  }

  unordered_map<bitset<nbit__>, size_t> blex;
  for (auto& ispace : *sdet->stringspaceb()) {
    for (auto& bbit : *ispace) blex[bbit] = ispace->lexical_zero(bbit);
  }

  auto lower_ras = [&sdet, &alex, &blex] (shared_ptr<const RASBlock<double>> sblock, shared_ptr<RASBlock<double>> tblock, const int nstart, const int nfence) {
    const size_t lb = sblock->lenb();
    double* odata = tblock->data();
    for (auto& abit : *tblock->stringsa()) {
      for (auto& bbit : *tblock->stringsb()) {
        for ( int i = nstart; i < nfence; ++i) {
          if (abit[i] || !bbit[i]) continue;
          bitset<nbit__> sabit = abit; sabit.set(i);
          bitset<nbit__> sbbit = bbit; sbbit.reset(i);

          const double phase = static_cast<double>(-1 * sdet->sign<0>(sabit, i) * sdet->sign<1>(sbbit,i));

          *odata += phase * sblock->element( blex[sbbit] + alex[sabit] * lb );
        }
        ++odata;
      }
    }
  };

  // The important thing to notice is that for all orbitals in a single RAS space, each block in the source is sent to a single block in target
  for (auto& iblock : out->blocks()) {
    if (!iblock) continue;
    const int nha = iblock->stringsa()->nholes();
    const int nhb = iblock->stringsb()->nholes();
    const int npa = iblock->stringsa()->nparticles();
    const int npb = iblock->stringsb()->nparticles();
    const int n2a = iblock->stringsa()->nele2();
    const int n2b = iblock->stringsb()->nele2();

    if ( (ras1 > 0) && (nhb < ras1) && (nha > 0) ) lower_ras(this->block(nha-1,nhb+1,npa,npb), iblock, 0, ras1);
    if ( (ras2 > 0) && (n2b > 0) && (n2a < ras2) ) lower_ras(this->block(nha, nhb, npa, npb), iblock, ras1, ras1 + ras2);
    if ( (ras3 > 0) && (npb > 0) && (npa < ras3) ) lower_ras(this->block(nha, nhb, npa+1, npb-1), iblock, ras1+ras2, ras1+ras2+ras3);
  }

  return out;
}

// S_+ = \sum_i i^dagger_alpha i_beta
template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert( (tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1) );
  auto out = make_shared<RASCivec>(tdet);

  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // maps bits to their local offsets
  unordered_map<bitset<nbit__>, size_t> alex;
  for (auto& ispace : *det_->stringspacea()) {
    for (auto& abit : *ispace) alex[abit] = ispace->lexical_zero(abit);
  }

  unordered_map<bitset<nbit__>, size_t> blex;
  for (auto& ispace : *det_->stringspaceb()) {
    for (auto& bbit : *ispace) blex[bbit] = ispace->lexical_zero(bbit);
  }

  auto raise_ras = [&sdet, &alex, &blex] (shared_ptr<const RASBlock<double>> sblock, shared_ptr<RASBlock<double>> tblock, const int nstart, const int nfence) {
    const size_t lb = sblock->lenb();
    double* odata = tblock->data();
    for (auto& abit : *tblock->stringsa()) {
      for (auto& bbit : *tblock->stringsb()) {
        for ( int i = nstart; i < nfence; ++i) {
          if (!abit[i] || bbit[i]) continue;
          bitset<nbit__> sabit = abit; sabit.reset(i);
          bitset<nbit__> sbbit = bbit; sbbit.set(i);

          const double phase = static_cast<double>(sdet->sign<0>(sabit, i) * sdet->sign<1>(sbbit,i));

          *odata += phase * sblock->element( blex[sbbit] + alex[sabit] * lb );
        }
        ++odata;
      }
    }
  };

  // The important thing to notice is that for all orbitals in a single RAS space, each block in the source is sent to a single block in target
  for (auto& iblock : out->blocks()) {
    if (!iblock) continue;
    const int nha = iblock->stringsa()->nholes();
    const int nhb = iblock->stringsb()->nholes();
    const int npa = iblock->stringsa()->nparticles();
    const int npb = iblock->stringsb()->nparticles();
    const int n2a = iblock->stringsa()->nele2();
    const int n2b = iblock->stringsb()->nele2();

    if ( (ras1 > 0) && (nha < ras1) && (nhb > 0) ) raise_ras(this->block(nha+1,nhb-1,npa,npb), iblock, 0, ras1);
    if ( (ras2 > 0) && (n2a > 0) && (n2b < ras2) ) raise_ras(this->block(nha, nhb, npa, npb), iblock, ras1, ras1 + ras2);
    if ( (ras3 > 0) && (npa > 0) && (npb < ras3) ) raise_ras(this->block(nha, nhb, npa-1, npb+1), iblock, ras1+ras2, ras1+ras2+ras3);
  }

  return out;
}

template<> void RASCivector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<RASCivec> S2 = spin();
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
