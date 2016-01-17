//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/denomtask.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef BAGEL_RAS_DENOMTASK_H
#define BAGEL_RAS_DENOMTASK_H

namespace bagel { namespace RAS {

/** Header for task to compute denominator for a RAS block. Compute function implemented in <src/ci/ras/rasci_denom.cc>a */
struct DenomTask {
  double* const data_;
  const std::bitset<nbit__> abit_;
  std::shared_ptr<const RASString> stringb_;
  const double* const jop_;
  const double* const kop_;
  const double* const h_;

  DenomTask(double* o, std::bitset<nbit__> ia, std::shared_ptr<const RASString> sb, double* j, double* k, double* h) :
    data_(o), abit_(ia), stringb_(sb), jop_(j), kop_(k), h_(h) {}

  void compute();
};

} }

#endif
