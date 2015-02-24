//
// BAGEL - Parallel electron correlation program.
// Filename: asd/asd_compute_rdm.hpp
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_COMPUTE_RDM_H
#define BAGEL_ASD_COMPUTE_RDM_H

template <class VecType>
void ASD<VecType>::compute_rdm12() {
  const int norbA = dimer_->active_refs().first->nact();
  const int norbB = dimer_->active_refs().second->nact();

  if (rdm1_av_ == nullptr && nstates_ > 1) {
    rdm1_av_ = std::make_shared<RDM<1>>(norbA+norbB);
    rdm2_av_ = std::make_shared<RDM<2>>(norbA+norbB);
  }
  if (nstates_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }

  compute_rdm12_dimer(); //allocation takes place

  compute_rdm12_monomer();

  //State-average RDM
  if (nstates_ != 1) {
    for (int i = 0; i != nstates_; ++i) {
      rdm1_av_->ax_plus_y(weight_[i], rdm1_[i]);
      rdm2_av_->ax_plus_y(weight_[i], rdm2_[i]);
    }  
  } else {                                      
    rdm1_av_ = rdm1_[0];
    rdm2_av_ = rdm2_[0];
  }
}


template <class VecType>
void ASD<VecType>::compute_rdm12_monomer() {
/*with in Monomer subspace
   ApBp
AB [ ]
*/
  std::cout << "monomer_rdm enetered" << std::endl;

  for (auto& AB : subspaces_) { //diagonal dimer subspace
    const int ioff = AB.offset();
    auto offset = std::make_pair(ioff,ioff);

    std::shared_ptr<const VecType> ccvecA = AB.template ci<0>(); //Dvec, RASDvec, DistDvec, DistRASDvec
    std::shared_ptr<const VecType> ccvecB = AB.template ci<1>();

    // A Ap B Bp
    std::array<VecType,4> fourvecs { ccvecA, ccvecA, ccvecB, ccvecB };
    std::array<MonomerKey,4> keys {AB.template monomerkey<0>(), AB.template monomerkey<0>(), 
                                   AB.template monomerkey<1>(), AB.template monomerkey<1>()};
    std::cout << keys[0].nstates() << std::endl; //TODO:delete

    for(int i = 0; i != nstates_; ++i) {
      std::shared_ptr<RDM<1>> rdm1;
      std::shared_ptr<RDM<2>> rdm2;
      std::tie(rdm1,rdm2) = compute_rdm12_monomer(i, offset, fourvecs);
      if (rdm1) *rdm1_[i] += *rdm1; 
      if (rdm2) *rdm2_[i] += *rdm2;
    }
  }
}

#endif

#endif
