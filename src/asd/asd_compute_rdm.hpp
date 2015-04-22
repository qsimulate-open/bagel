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
  Timer rdmtime;
  std::cout << std::endl << " ===== ASD RDM Computation ==== " << std::endl;
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
  std::cout << "  o Dimer RDM - " << std::setw(9) << std::fixed << std::setprecision(2) << rdmtime.tick() << std::endl;

#if 1 //Monomer
  //1 monomer RDM is calculated with standard algorithm
  //0 calculated with Gamma tensors (then make sure Monomer switch in other files are set to 1)
  compute_rdm12_monomer();
  std::cout << "  o Momer RDM - " << std::setw(9) << std::fixed << std::setprecision(2) << rdmtime.tick() << std::endl;
#endif

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

  //debug
  for (int i = 0; i != nstates_; ++i) {
    debug_rdm(rdm1_[i], rdm2_[i], i, /*mute*/false);
    debug_energy(rdm1_[i], rdm2_[i], i, /*mute*/false);
  }
  std::cout << "  o Debug RDM - " << std::setw(9) << std::fixed << std::setprecision(2) << rdmtime.tick() << std::endl;
  std::cout << " ============================== " << std::endl;
}


template <class VecType>
void ASD<VecType>::compute_rdm12_monomer() {
  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();

  std::vector<std::shared_ptr<RDM<1>>> rdm1A; rdm1A.resize(nstates_);
  std::vector<std::shared_ptr<RDM<2>>> rdm2A; rdm2A.resize(nstates_);
  std::vector<std::shared_ptr<RDM<1>>> rdm1B; rdm1B.resize(nstates_);
  std::vector<std::shared_ptr<RDM<2>>> rdm2B; rdm2B.resize(nstates_);

  for (auto& AB : subspaces_) { //diagonal dimer subspace
    const int ioff = AB.offset();
    std::shared_ptr<const VecType> A = AB.template ci<0>();
    std::shared_ptr<const VecType> B = AB.template ci<1>();
    const int nstA = A->ij(); //Dvec size
    const int nstB = B->ij();
    assert(nstA == AB.nstatesA());
    assert(nstB == AB.nstatesB());

    //MonomerA i.e. delta_JJ'
    for (int i = 0; i != nstA; ++i) { // <I|
      for (int ip = 0; ip != nstA; ++ip) { // |I'>

        std::shared_ptr<RDM<1>> rdm1;
        std::shared_ptr<RDM<2>> rdm2;
        std::tie(rdm1,rdm2) = compute_rdm12_monomer(A, i, ip); // <I'|E(op)|I>

        for (int kst = 0; kst != nstates_; ++kst) {

          double csum = 0.0; //coeff sum
          for (int j = 0; j != nstB; ++j) { // delta_J'J
            const int ij  = i  + (j*nstA);
            const int ijp = ip + (j*nstA);
            assert(ij < adiabats_->ndim());
            assert(ijp < adiabats_->ndim());
            csum += adiabats_->element(ioff+ij,kst) * adiabats_->element(ioff+ijp,kst); //cf. dimerindex()
          }

          if (!rdm1A[kst]) rdm1A[kst] = std::make_shared<RDM<1>>(nactA);
          if (!rdm2A[kst]) rdm2A[kst] = std::make_shared<RDM<2>>(nactA);
          rdm1A[kst]->ax_plus_y(csum, rdm1);
          rdm2A[kst]->ax_plus_y(csum, rdm2);

        }
      }
    }

    //MonomerB i.e. delta_II'
    for (int j = 0; j != nstB; ++j) { //<J|
      for (int jp = 0; jp != nstB; ++jp) { // |J'>

        std::shared_ptr<RDM<1>> rdm1;
        std::shared_ptr<RDM<2>> rdm2;
        std::tie(rdm1,rdm2) = compute_rdm12_monomer(B, j, jp); // <J'|E(op)|J>

        for (int kst = 0; kst != nstates_; ++kst) {

          double csum = 0.0;
          for (int i = 0; i != nstA; ++i) { // delta_I'I
            const int ij  = i + (j*nstA);
            const int ijp = i + (jp*nstA);
            assert(ij < adiabats_->ndim());
            assert(ijp < adiabats_->ndim());
            csum += adiabats_->element(ioff+ij,kst) * adiabats_->element(ioff+ijp,kst);
          }

          if (!rdm1B[kst]) rdm1B[kst] = std::make_shared<RDM<1>>(nactB);
          if (!rdm2B[kst]) rdm2B[kst] = std::make_shared<RDM<2>>(nactB);
          rdm1B[kst]->ax_plus_y(csum, rdm1);
          rdm2B[kst]->ax_plus_y(csum, rdm2);

        }
      }
    } //B
  } //subspaces

  for (int istate = 0; istate != nstates_; ++istate) {
    auto rdm1 = std::make_shared<RDM<1>>(nactA+nactB);
    {//A
      auto low = {0,0};
      auto up  = {nactA,nactA};
      auto outv = make_rwview(rdm1->range().slice(low,up), rdm1->storage());
      copy(rdm1A[istate]->begin(), rdm1A[istate]->end(), outv.begin());
    }
    {//B
      auto low = {nactA,nactA};
      auto up  = {nactA+nactB,nactA+nactB};
      auto outv = make_rwview(rdm1->range().slice(low,up), rdm1->storage());
      copy(rdm1B[istate]->begin(), rdm1B[istate]->end(), outv.begin());
    }
    auto rdm2 = std::make_shared<RDM<2>>(nactA+nactB);
    {//A
      auto low = {0,0,0,0};
      auto up  = {nactA,nactA,nactA,nactA};
      auto outv = make_rwview(rdm2->range().slice(low,up), rdm2->storage());
      copy(rdm2A[istate]->begin(), rdm2A[istate]->end(), outv.begin());
    }
    {//B
      auto low = {nactA,nactA,nactA,nactA};
      auto up  = {nactA+nactB,nactA+nactB,nactA+nactB,nactA+nactB};
      auto outv = make_rwview(rdm2->range().slice(low,up), rdm2->storage());
      copy(rdm2B[istate]->begin(), rdm2B[istate]->end(), outv.begin());
    }

    //update rdms
    *rdm1_[istate] += *rdm1;
    *rdm2_[istate] += *rdm2;
  }

}

#endif

#endif
