//
// BAGEL - Parallel electron correlation program.
// Filename: asd/asd_rdm.cc
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

#include <src/asd/asd.h>
#include <src/asd/state_tensor.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace btas;

template <class VecType>
void ASD<VecType>::compute_rdm12() {
  const int norbA = dimer_->active_refs().first->nact();
  const int norbB = dimer_->active_refs().second->nact();

  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norbA+norbB);
    rdm2_av_ = make_shared<RDM<2>>(norbA+norbB);
  }
  if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }

  // compute transformed gammas (J',J,zeta)
  StateTensor st(adiabats_, subspaces_base());
  st.print();

  for (int i = 0; i != nstate_; ++i) compute_rdm12(i);
  

}


template <class VecType>
void ASD<VecType>::compute_rdm12(const int istate) {
  const int norbA = dimer_->active_refs().first->nact();
  const int norbB = dimer_->active_refs().second->nact();

  auto rdm1 = make_shared<RDM<1>>(norbA+norbB);
  auto rdm2 = make_shared<RDM<2>>(norbA+norbB);

  += compute_rdm12_monomer(istate);

  += compute_rdm12_dimer(istate);

  symmetrize_RDM(); //TODO: apply to each RDM
  debug_RDM(); 
  debug_energy();
  
  // setting to private members.
  rdm1_[istate] = rdm1;
  rdm2_[istate] = rdm2;
  if (nstate_ != 1) {
    rdm1_av_->ax_plus_y(weight_[istate], rdm1);
    rdm2_av_->ax_plus_y(weight_[istate], rdm2);
  } else {
    rdm1_av_ = rdm1;
    rdm2_av_ = rdm2;
  }

}


template <class VecType>
std::tuple<std::shared_ptr<RDM<1>>,std::shared_ptr<RDM<2>>> ASD<VecType>::compute_rdm12_monomer(const int istate) {
/* with in Monomer subspace
   ApBp
AB [ ]
TODO/ This function is written in such way to facilitate extension to offdiagonal subspaces
*/
  std::cout << "monomer_rdm enetered" << std::endl;
  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();

//int istate = 0; // ground state
  //RDM
  onerdm_ = std::make_shared<RDM<1>>(nactA+nactB);
  twordm_ = std::make_shared<RDM<2>>(nactA+nactB);
  

  // Diagonal Dimer subspaces
  int isub = 0;
  for (auto& AB : subspaces_) {
    //prints the CI vectors of each dimerspace
    isub++;
    const int ioff = AB.offset();
    std::cout << "Dimer Subspace #" << isub << ": offset = " << ioff << std::endl;
    auto offset = std::make_pair(ioff,ioff);

    std::shared_ptr<const VecType> ccvecA = AB.template ci<0>(); //Dvec, RASDvec, DistDvec, DistRASDvec
    std::cout << "CI vectors of Monomer A" << std::endl;
    ccvecA->print(); //in case of Dvec (fci/dvec.h)
    std::shared_ptr<const VecType> ccvecB = AB.template ci<1>();
    std::cout << "CI vectors of Monomer B" << std::endl;
    ccvecB->print();

    std::cout << "Dimension of subspace = " << ccvecA->ij() << " x " << ccvecB->ij() << std::endl;

    std::shared_ptr<RDM<1>> r1;
    std::shared_ptr<RDM<2>> r2;

    // A Ap B Bp
    std::array<VecType,4> fourvecs { ccvecA, ccvecA, ccvecB, ccvecB };
    std::array<MonomerKey,4> keys {AB.template monomerkey<0>(), AB.template monomerkey<0>(), 
                                   AB.template monomerkey<1>(), AB.template monomerkey<1>()};
    std::cout << keys[0].nstates() << std::endl; //TODO:delete

    tie(r1,r2) = compute_rdm12_monomer(offset, fourvecs); /*implemented in cas, ras, distcas, distras*/
    if(r1) *onerdm_ += *r1;
    if(r2) *twordm_ += *r2;
  }
}



//SKIP RDM34
#if 0
  //3&4RDM
  threerdm_ = std::make_shared<RDM<3>>(nactA+nactB);
  threerdm_->zero();
  initialize_4RDM();
    //3&4RDM
    std::shared_ptr<RDM<3>> r3;
    std::shared_ptr<RDM<4>> r4A;
    std::shared_ptr<RDM<4>> r4B;
    tie(r3,r4A,r4B) = compute_rdm34_monomer(offset, fourvecs);
    if(r3) *threerdm_ += *r3;
    if(r4A) {
      btas::CRange<2> range(nactA*nactA*nactA*nactA*nactA*nactA*nactA*nactA,1);
      const MatView matv(btas::make_view(range, r4A->storage()), /*localized*/true);
      auto mat = std::make_shared<Matrix>( matv );
      *fourrdmparts_.at("monomerA") += *mat;
    }
    if(r4B) {
      btas::CRange<2> range(nactB*nactB*nactB*nactB*nactB*nactB*nactB*nactB,1);
      const MatView matv(btas::make_view(range, r4B->storage()), /*localized*/true);
      auto mat = std::make_shared<Matrix>( matv );
      *fourrdmparts_.at("monomerB") += *mat;
    }
  }
 
/*
  // Offdiagonal Dimer subspaces
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    for (auto jAB = subspaces_.begin(); jAB != subspaces_.end(); ++jAB) {
      if (iAB == jAB) continue;
      //TODO Lower-triangular (i<->j)
      auto offset = std::make_pair(iAB->offset(),jAB->offset());

      std::shared_ptr<const VecType> ccvecA  = iAB->template ci<0>(); // <I'(A)|//Dvec, RASDvec, DistDvec, DistRASDvec
      std::shared_ptr<const VecType> ccvecB  = iAB->template ci<1>(); // <J'(B)|
      std::shared_ptr<const VecType> ccvecAp = jAB->template ci<0>(); // |I(A)> 
      std::shared_ptr<const VecType> ccvecBp = jAB->template ci<1>(); // |J(B)>
      // A Ap B Bp
      std::array<VecType,4> fourvecs { ccvecA, ccvecAp, ccvecB, ccvecBp };

      std::shared_ptr<RDM<1>> r1;
      std::shared_ptr<RDM<2>> r2;
      Coupling term_type = coupling_type(*iAB, *jAB);
      switch(term_type) {
        case Coupling::diagonal :
          std::cout << "Offdiagonal subspace, diagonal coupling, offset = " << std::get<0>(offset) << " and " << std::get<1>(offset) << std::endl;
          tie(r1,r2) = compute_rdm12_monomer(offset, fourvecs);
          break;
        default : //do nothing
          break;
      }
      if (r1) {
        *onerdm_ += *r1; //TODO : *2 if (lower triangular is used)
      }
      if (r2) {
        *twordm_ += *r2;
      }

    }
  }
*/

  //PRINT
  std::cout << "!@# Monomer RDM print" << std::endl;
  std::cout << "Active space: A(" << nactA << "), B(" << nactB << ")" << std::endl;
  onerdm_->print(1.0e-6);

  //APPROX 2RDM
  approx2rdm_ = std::make_shared<RDM<2>>(*twordm_);
}

#endif
