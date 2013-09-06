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


#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include <bitset>

#include <src/ras/civector.h>
#include <src/parallel/resources.h> // This is maybe only necessary because I am missing something else?
#include <src/util/constants.h>
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

      SpinTask(const std::bitset<nbit__> t, const RASCivector<double>* th, shared_ptr<RASCivector<double>> o, shared_ptr<const RASDeterminants> d) :
        target_(t), this_(th), out_(o), det_(d) {}

      void compute() {
        const int norb = det_->norb();

        unique_ptr<double[]> source(new double[det_->lenb()]);

        for (auto& iter : det_->phia(det_->lexical<0>(target_))) {
          const int ii = iter.ij / norb;
          const int jj = iter.ij % norb;
          bitset<nbit__> mask1; mask1.set(ii); mask1.set(jj);
          bitset<nbit__> mask2; mask2.set(ii);

          bitset<nbit__> maskij; maskij.set(ii); maskij.flip(jj);

          fill_n(source.get(), det_->lenb(), 0.0);
          vector<shared_ptr<RASBlock<double>>> sourceblocks = out_->allowed_blocks<0>(det_->stringa(iter.source));
          for (auto& iblock : sourceblocks) {
            const int offset = iblock->stringb()->offset();
            copy_n(&this_->element(iblock->stringb()->strings(0), det_->stringa(iter.source)), iblock->lenb(), source.get()+offset);
          }

          vector<shared_ptr<RASBlock<double>>> allowedblocks = out_->allowed_blocks<0>(target_);
          for (auto& iblock : allowedblocks) {
            double* outelement = &out_->element(iblock->stringb()->strings(0), target_);
            for (auto& btstring : *iblock->stringb()) {
              if ( ((btstring & mask1) ^ mask2).none() ) { // equivalent to "btstring[ii] && (ii == jj || !btstring[jj])"
                const bitset<nbit__> bsostring = btstring ^ maskij;
                if (det_->allowed(det_->stringa(iter.source), bsostring))
                  *outelement -= static_cast<double>(iter.sign * det_->sign(bsostring, ii, jj)) * source[det_->lexical<1>(bsostring)];
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

  TaskQueue<RAS::SpinTask> tasks(det_->stringa().size());

  for (auto& istring : det_->stringa()) {
    tasks.emplace_back(istring, this, out, det_);
  }

  tasks.compute();

  const double sz = static_cast<double>(det_->nspin()) * 0.5;
  const double fac = sz*sz + sz + static_cast<double>(det_->neleb());

  out->ax_plus_y(fac, *this);

  return out;
}

// S_-
template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(shared_ptr<const RASDeterminants>) const {
  assert(false);
  return shared_ptr<RASCivector<double>>();
}

// S_+
template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(shared_ptr<const RASDeterminants>) const {
  assert(false);
  return shared_ptr<RASCivector<double>>();
}

template<> void RASCivector<double>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<RASCivec> S2 = spin();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while( fabs(actual_expectation - pure_expectation) > thresh ) {
    cout << "<S^2> = " << actual_expectation << endl;
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
