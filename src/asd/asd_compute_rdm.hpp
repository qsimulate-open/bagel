//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/asd_compute_rdm.hpp
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

  compute_rdm12_monomer();
  std::cout << "  o Monomer RDM - " << std::setw(9) << std::fixed << std::setprecision(2) << rdmtime.tick() << std::endl;

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

  //print information
  if (print_info_) {
    for (int i = 0; i != nstates_; ++i) {
      print_rdm_info(rdm1_[i], rdm2_[i], i);
      print_energy_info(rdm1_[i], rdm2_[i], i);
    }
  }

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

  //allocate first (filled with zero)
  for (int kst = 0; kst != nstates_; ++kst) {
    rdm1A[kst] = std::make_shared<RDM<1>>(nactA);
    rdm2A[kst] = std::make_shared<RDM<2>>(nactA);
    rdm1B[kst] = std::make_shared<RDM<1>>(nactB);
    rdm2B[kst] = std::make_shared<RDM<2>>(nactB);
  }

  const int n = subspaces_.size() * nstates_;
  const int numproc = mpi__->size();
  const int myid = mpi__->rank();

  std::vector<std::pair<int,int>> proc;
  for (int i = 0; i != subspaces_.size(); ++i)
    for (int j = 0; j != nstates_; ++j)
      proc.emplace_back(i,j);
  assert(proc.size() == n);

  int mystart, myend;
  if (n >= numproc) {
    mystart = (n / numproc) * myid + ((n % numproc) < myid ? (n % numproc) : myid);
    myend = mystart + (n / numproc) + ((n % numproc) > myid);
  }
  else {
    mystart = myid < n ? myid : -1;
    myend = myid < n ? (myid + 1) : 0;
  }

  //diagonal dimer subspace
  for (int job = mystart; job != myend; ++job) {

    if (job < 0) break; //MPI workers without allocated jobs, fill zero and wait for reduction

    const int isub = proc[job].first;
    const int kst = proc[job].second;

    auto& AB = subspaces_[isub];
    std::shared_ptr<const VecType> A = AB.template ci<0>(); //CASDvec, RASDvec
    std::shared_ptr<const VecType> B = AB.template ci<1>();
    const int ioff = AB.offset();
    const int nstA = A->ij(); //XDvec size
    const int nstB = B->ij();
    assert(nstA == AB.nstatesA());
    assert(nstB == AB.nstatesB());

    std::shared_ptr<RDM<1>> rdm1;
    std::shared_ptr<RDM<2>> rdm2;

    {//MonomerA i.e. delta_JJ'
      auto dket = contract_I(A, adiabats_, ioff, nstA, nstB, kst);
      for (int j = 0; j != nstB; ++j) {
        std::tie(rdm1, rdm2) = compute_rdm12_monomer(dket, j);
        rdm1A[kst]->ax_plus_y(1.0, rdm1);
        rdm2A[kst]->ax_plus_y(1.0, rdm2);
      }
    }

    {//MonomerB i.e. delta_II'
      auto dket = contract_J(B, adiabats_, ioff, nstA, nstB, kst);
      for (int i = 0; i != nstA; ++i) {
        std::tie(rdm1, rdm2) = compute_rdm12_monomer(dket, i);
        rdm1B[kst]->ax_plus_y(1.0, rdm1);
        rdm2B[kst]->ax_plus_y(1.0, rdm2);
      }
    }
  } //subspaces

#ifdef HAVE_MPI_H
  for (int kst = 0; kst != nstates_; ++kst) {
    mpi__->allreduce(rdm1A[kst]->data(), rdm1A[kst]->size());
    mpi__->allreduce(rdm2A[kst]->data(), rdm2A[kst]->size());
    mpi__->allreduce(rdm1B[kst]->data(), rdm1B[kst]->size());
    mpi__->allreduce(rdm2B[kst]->data(), rdm2B[kst]->size());
  }
#endif

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
