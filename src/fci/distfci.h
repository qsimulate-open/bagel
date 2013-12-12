//
// BAGEL - Parallel electron correlation program.
// Filename: distfci.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __BAGEL_FCI_DISTFCI_H
#define __BAGEL_FCI_DISTFCI_H

#include <memory>
#include <vector>
#include <src/fci/fci.h>
#include <src/fci/space.h>

namespace bagel {

// Parallel FCI based on Harrison-Zarrabian algorithm.
// The algorithm is basically the same as written in
// Z. Gan and R. J. Harrison, SC '05: Proceedings of the 2005 ACM/IEEE conference on Supercomputing
//
// The implementation is based on the HarrisonZarrabian class written by Shane Parker.

class DistFCI : public Method {
  protected:
    std::shared_ptr<Space> space_;
    std::shared_ptr<DistCivec> denom_;

    // Options
    int max_iter_;
    int davidsonceiling_;
    double thresh_;
    double print_thresh_;

    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;

    int nstate_;

    // extra
    std::shared_ptr<const Determinants> det_;

    // results
    std::vector<double> energy_;
    std::shared_ptr<DistDvec> cc_;
    std::shared_ptr<MOFile> jop_;

    void common_init();
    void print_header() const;

    // const_denom function here only makes a denom for local data of DistCivec.
    void const_denom();

    // denominator
    void generate_guess(const int nspin, const int nstate, std::vector<std::shared_ptr<DistCivec>> out);

    // Determinant seeds in parallel
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet);

  public:
    // this constructor is ugly... to be fixed some day...
    DistFCI(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
            const int ncore = -1, const int nocc = -1, const int nstate = -1);

    // FCI compute function using DistCivec
    void compute() override;

    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    double core_energy() const { return jop_->core_energy(); }

    int nij() const { return (norb_*(norb_+1))/2; }

    void update(std::shared_ptr<const Coeff>);

    std::shared_ptr<const Determinants> det() const { return det_; }
    std::shared_ptr<const MOFile> jop() const { return jop_; }
    std::shared_ptr<const DistCivec> denom() const { return denom_; }
    std::shared_ptr<const DistDvec> civectors() const { return cc_; }

    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    std::shared_ptr<const Reference> conv_to_ref() const override { return std::shared_ptr<const Reference>(); }
};

}

#endif

