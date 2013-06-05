
// BAGEL - Parallel electron correlation program.
// Filename: dimer_jop.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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



#ifndef __BAGEL_DIMER_JOP_H
#define __BAGEL_DIMER_JOP_H

#include <src/fci/mofile.h>

namespace bagel {

class DimerJop : public Jop {
  protected:
    std::pair<std::unique_ptr<double[]>, std::unique_ptr<double[]>> monomer_mo1es_;
    std::pair<std::unique_ptr<double[]>, std::unique_ptr<double[]>> monomer_mo2es_;

    std::shared_ptr<Matrix> cross_mo1e_;

  public:
    DimerJop(const std::shared_ptr<const Reference> ref, const int nstart, const int nfenceA, const int nfenceB,
      std::shared_ptr<const Coeff> coeff); // note that in DimerJop, I'm forcing a HZ Jop
    ~DimerJop() {};

    // Functions to kind of make DimerJop behave like a pair of shared_ptr<Jop>
    double* mo1e_first() const { return monomer_mo1es_.first.get(); };
    double* mo1e_second() const { return monomer_mo1es_.second.get(); };

    double* mo2e_first() const { return monomer_mo2es_.first.get(); };
    double* mo2e_second() const { return monomer_mo2es_.second.get(); };

    std::shared_ptr<Matrix> cross_mo1e() const { return cross_mo1e_; }
};

}

#endif
