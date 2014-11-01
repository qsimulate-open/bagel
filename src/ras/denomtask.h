//
// BAGEL - Parallel electron correlation program.
// Filename: ras/denomtask.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef BAGEL_RAS_DENOMTASK_H
#define BAGEL_RAS_DENOMTASK_H

namespace bagel { namespace RAS {

/** Header for task to compute denominator for a RAS block. Compute function implemented in <src/ras/rasci_denom.cc>a */
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
