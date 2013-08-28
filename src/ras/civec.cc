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

        for (auto& iter : det_->phia(det_->lexical<0>(target_))) {
          const int jj = iter.ij / norb;
          const int ii = iter.ij - norb * jj;

          vector<RASBlock<double>> allowedblocks = out_->allowed_blocks<0>(target_);
          for (auto& iblock : allowedblocks) {
            for (auto& btstring : *iblock.stringb()) {
              if ( !btstring[jj] || ( ii != jj && btstring[ii] ) ) continue;
              bitset<nbit__> bsostring = btstring; bsostring.reset(jj); bsostring.set(ii);
              out_->element(target_, btstring) -= static_cast<double>(iter.sign * det_->sign(bsostring, ii, jj)) * this_->element(iter.source, bsostring);
            }
          }
        }
      }
    };
  }
}

// Returns S^2 | civec >
// The important portion of S^2 is S_-S_+ = nbeta - \sum_{ij} j_alpha^dagger i_alpha i_beta^dagger j_alpha
template<>
shared_ptr<RASCivector<double>> RASCivector<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);

  vector<RAS::SpinTask> tasks;
  tasks.reserve( det_->stringa().size() );

  for (auto& istring : det_->stringa())
    tasks.emplace_back(istring, this, out, det_);

  TaskQueue<RAS::SpinTask> tq(tasks);
  tq.compute(resources__->max_num_threads());

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

  const double expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<RASCivec> S2 = spin();

  int k = nspin + 2;
  while( fabs(dot_product(*S2) - expectation) > thresh ) {
    if ( k > max_spin ) throw runtime_error("Spin decontamination failed.");

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();

    k += 2;
  }
}
