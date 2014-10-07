#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_MONOMER_RDM_H
#define BAGEL_ASD_MONOMER_RDM_H


#include <src/asd/asd_base.h>

//#include <src/asd/coupling.h>

//***************************************************************************************************************
template <class VecType>
void 
ASD<VecType>::monomer_rdm () {
/*
@called from:
ASD::asd_compute()
@note:
with in Monomer subspace
   ApBp
AB [ ]
TODO/ This function is written in such way to facilitate extension to offdiagonal subspaces
*/
//***************************************************************************************************************
  std::cout << "monomer_rdm enetered" << std::endl;
  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();

//int istate = 0; // ground state
  //RDM
  onerdm_ = std::make_shared<RDM<1>>(nactA+nactB);
  twordm_ = std::make_shared<RDM<2>>(nactA+nactB);
  onerdm_->zero();
  twordm_->zero();

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

    tie(r1,r2) = compute_rdm12_monomer(offset, fourvecs);
    if(r1) *onerdm_ += *r1;
    if(r2) *twordm_ += *r2;
  }

/*
  // Offdiagonal Dimer subspaces
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      //Lower-triangular (i<->j)
      auto offset = std::make_pair(jAB->offset(),iAB->offset());

      std::shared_ptr<const VecType> ccvecA  = jAB->template ci<0>(); //Dvec, RASDvec, DistDvec, DistRASDvec
      std::shared_ptr<const VecType> ccvecAp = iAB->template ci<0>();
      std::shared_ptr<const VecType> ccvecB  = jAB->template ci<1>();
      std::shared_ptr<const VecType> ccvecBp = iAB->template ci<1>();
      // A Ap B Bp
      std::array<VecType,4> fourvecs { ccvecA, ccvecAp, ccvecB, ccvecBp };

      std::shared_ptr<RDM<1>> r1;
      std::shared_ptr<RDM<2>> r2;
      Coupling term_type = coupling_type(*jAB, *iAB);
      switch(term_type) {
        case Coupling::diagonal :
          std::cout << "Offdiagonal subspace, diagonal coupling, offset = " << std::get<0>(offset) << " and " << std::get<1>(offset) << std::endl;
          tie(r1,r2) = compute_rdm12_monomer(offset, fourvecs);
          break;
        default : //do nothing
          break;
      }
      if (r1) {
        *onerdm_ += *r1; //TODO : *2
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

}

#endif

#endif

