//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/form_sigma.cc
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

#include <map>

#include <src/asd_dmrg/form_sigma.h>
#include <src/ras/form_sigma.h>
#include <src/ras/apply_operator.h>

using namespace std;
using namespace bagel;

vector<shared_ptr<ProductRASCivec>> FormSigmaProdRAS::operator()(vector<shared_ptr<ProductRASCivec>>& ccvec,
                     shared_ptr<const BlockOperators> blockops, shared_ptr<const DimerJop> jop, const vector<bool>& conv) const {

  const int nstate = ccvec.size();
  vector<shared_ptr<ProductRASCivec>> sigmavec;
  for_each(ccvec.begin(), ccvec.end(), [&sigmavec] (shared_ptr<const ProductRASCivec> c) { sigmavec.push_back(c->clone()); });

  for (int istate = 0; istate != nstate; ++istate) {
#ifdef NOPENOPE_HAVE_MPI_H
    if ( istate % mpi__->size() == mpi__->rank() ) {
#endif
      Timer pdebug;
      shared_ptr<const ProductRASCivec> cc = ccvec.at(istate);
      shared_ptr<ProductRASCivec> sigma = sigmavec.at(istate);

      // pure terms
      pure_block_and_ras(cc, sigma, blockops, jop);
      pdebug.tick_print("pure");

      interaction_terms(cc, sigma, blockops, jop);
      pdebug.tick_print("interaction");
#ifdef HAVE_MPI_H
    }
#endif
  }

#ifdef NOPENOPE_HAVE_MPI_H
  for (int istate = 0; istate != nstate; ++istate) {
    if (!conv[istate])
      sigmavec.at(istate)->broadcast(istate & mpi__->size()); // This broadcast function is not yet written
  }
#endif

  return sigmavec;
}

void FormSigmaProdRAS::pure_block_and_ras(shared_ptr<const ProductRASCivec> cc, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops, shared_ptr<const DimerJop> jop) const {
  Timer ptime;
  for (auto& sector : sigma->sectors()) {
    // first prepare pure block part which will be a nsecstates x nsecstates matrix
    const int nsecstates = sector.second->nstates();
    const Matrix& pure_block = *blockops->ham(sector.first);
    assert(pure_block.ndim()==nsecstates && pure_block.mdim()==nsecstates);

    shared_ptr<RASBlockVectors> sigma_sector = sector.second;
    shared_ptr<const RASBlockVectors> cc_sector = cc->sector(sector.first);
    dgemm_("N","T", sigma_sector->ndim(), sigma_sector->mdim(), sigma_sector->mdim(), 1.0, cc_sector->data(), cc_sector->ndim(), pure_block.data(), pure_block.ndim(),
                                                                                      0.0, sigma_sector->data(), sigma_sector->ndim());
    ptime.tick_print("pure_block");

    // now do individual form_sigmas for the RAS parts
    FormSigmaRAS form_pure_ras(batchsize_);
    for(int ist = 0; ist < nsecstates; ++ist)
      form_pure_ras(cc_sector->civec(ist), sigma_sector->civec(ist), jop->monomer_jop<0>());
    ptime.tick_print("pure_ras");
  }
}


void FormSigmaProdRAS::interaction_terms(shared_ptr<const ProductRASCivec> cc, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops, shared_ptr<const DimerJop> jop) const {
  Timer ptime;

  for (auto& isec : cc->sectors()) {
    shared_ptr<const RASBlockVectors> cc_sector = isec.second;
    const BlockKey cc_key = isec.first;

    // really ugly but verbose
    const bool do_aET    = cc->left()->contains(BlockKey(cc_key.nelea-1, cc_key.neleb  ));
    const bool do_bET    = cc->left()->contains(BlockKey(cc_key.nelea  , cc_key.neleb-1));
    const bool do_aHT    = cc->left()->contains(BlockKey(cc_key.nelea+1, cc_key.neleb  ));
    const bool do_bHT    = cc->left()->contains(BlockKey(cc_key.nelea  , cc_key.neleb+1));
    const bool do_aaET   = cc->left()->contains(BlockKey(cc_key.nelea-2, cc_key.neleb  ));
    const bool do_bbET   = cc->left()->contains(BlockKey(cc_key.nelea  , cc_key.neleb-2));
    const bool do_abET   = cc->left()->contains(BlockKey(cc_key.nelea-1, cc_key.neleb-1));
    const bool do_aaHT   = cc->left()->contains(BlockKey(cc_key.nelea+2, cc_key.neleb  ));
    const bool do_bbHT   = cc->left()->contains(BlockKey(cc_key.nelea  , cc_key.neleb+2));
    const bool do_abHT   = cc->left()->contains(BlockKey(cc_key.nelea+1, cc_key.neleb+1));
    const bool do_flipup = cc->left()->contains(BlockKey(cc_key.nelea-1, cc_key.neleb+1));
    const bool do_flipdn = cc->left()->contains(BlockKey(cc_key.nelea+1, cc_key.neleb-1));

    if (do_aET || do_aaET) {
      branch_1(cc_sector, sigma, blockops);
      ptime.tick_print("branch-1");
    }

    if (do_bET || do_bbET || do_abET) {
      branch_2(cc_sector, sigma, blockops);
      ptime.tick_print("branch-2");
    }

    if (do_aHT || do_aaHT || do_bHT) {
      branch_3(cc_sector, sigma, blockops);
      ptime.tick_print("branch-3");
    }

    if (do_bHT || do_aHT || do_bbHT || do_abHT) {
      branch_4(cc_sector, sigma, blockops);
      ptime.tick_print("branch-4");
    }

    // always compute branch 5 and 6
    branch_5(cc_sector, sigma, blockops);
    branch_6(cc_sector, sigma, blockops);
    ptime.tick_print("branch-5&6");

    if (do_flipup || do_aET) {
      branch_7(cc_sector, sigma, blockops);
      ptime.tick_print("branch-7");
    }

    if (do_flipdn || do_bET) {
      branch_8(cc_sector, sigma, blockops);
      ptime.tick_print("branch-8");
    }
  }
}

// should only enter this code if at least one of the terms will be computed
void FormSigmaProdRAS::branch_1(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();
  // S_alpha^+
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey singleETkey(cckey.nelea-1, cckey.neleb);
  const BlockKey doubleETkey(cckey.nelea-2, cckey.neleb);

  const bool do_single = sigma->left()->contains(singleETkey);
  const bool do_double = sigma->left()->contains(doubleETkey);
  // not sure how you could do a double but not a single, but that's a different problem
  assert(do_single || do_double);

  shared_ptr<const btas::Tensor3<double>> S = do_single ? blockops->S_a(singleETkey) : nullptr;
  shared_ptr<const btas::Tensor3<double>> P = do_double ? blockops->P_aa(cckey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> single_sector = do_single ? sigma->sector(singleETkey) : nullptr;
  shared_ptr<const RASDeterminants> single_det = do_single ? single_sector->det() : sigma->space()->det(singleETkey.nelea, singleETkey.neleb);

  const BlockInfo singlestate(singleETkey.nelea, singleETkey.neleb, nccstates);
  RASBlockVectors sector_r(single_det, singlestate);

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::CreateAlpha}, {r});
    if (do_single)
      dgemm_("N", "N", single_sector->ndim(), single_sector->mdim(), sector_r.mdim(), 1.0, sector_r.data(), sector_r.ndim(),
                                                          &(*S)(0,0,r), S->extent(0), 1.0, single_sector->data(), single_sector->ndim());
#if 0
    if (do_double) {
      // some extra stuff
    }
#endif
  }
}

// should only enter this code if at least one of the terms will be computed
void FormSigmaProdRAS::branch_2(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  // S_beta^+
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey bETkey(cckey.nelea, cckey.neleb-1);
  const BlockKey bbETkey(cckey.nelea, cckey.neleb-2);
  const BlockKey abETkey(cckey.nelea-1, cckey.neleb-1);

  const bool do_b = sigma->left()->contains(bETkey);
  const bool do_bb = sigma->left()->contains(bbETkey);
  const bool do_ab = sigma->left()->contains(abETkey);
  assert(do_b || do_bb || do_ab);

  shared_ptr<const btas::Tensor3<double>> S = do_b ? blockops->S_b(bETkey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_bb = do_bb ? blockops->P_bb(cckey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_ab = do_ab ? blockops->P_ab(cckey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> b_sector = do_b ? sigma->sector(bETkey) : nullptr;
  shared_ptr<const RASDeterminants> b_det = do_b ? b_sector->det() : sigma->space()->det(bETkey.nelea, bETkey.neleb);

  const BlockInfo bstate(bETkey.nelea, bETkey.neleb, nccstates);
  RASBlockVectors sector_r(b_det, bstate);

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::CreateBeta}, {r});
    if (do_b)
      dgemm_("N", "N", b_sector->ndim(), b_sector->mdim(), sector_r.mdim(), 1.0, sector_r.data(), sector_r.ndim(),
                                                &(*S)(0,0,r), S->extent(0), 1.0, b_sector->data(), b_sector->ndim());
#if 0
    if (do_bb) {
      // some extra stuff
    }

    if (do_ab) {
      // some extra stuff
    }
#endif
  }
}

// should only enter this code if at least one of the terms will be computed
void FormSigmaProdRAS::branch_3(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  // S_alpha
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey aHTkey(cckey.nelea+1, cckey.neleb);
  const BlockKey aaHTkey(cckey.nelea+2, cckey.neleb);
  const BlockKey bHTkey(cckey.nelea, cckey.neleb+1);

  const bool do_aHT = sigma->left()->contains(aHTkey);
  const bool do_aaHT = sigma->left()->contains(aaHTkey);
  const bool do_bHT = sigma->left()->contains(bHTkey);
  // not sure how you could do a aa but not a a, but that's a different problem
  assert(do_aHT || do_aaHT || do_bHT);

  shared_ptr<const btas::Tensor3<double>> S = do_aHT ? blockops->S_a(cckey) : nullptr;
  shared_ptr<const btas::TensorN<double,5>> Da = do_aHT ? blockops->D_a(cckey) : nullptr;
  shared_ptr<const btas::TensorN<double,5>> Db = do_bHT ? blockops->D_b(cckey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P = do_aaHT ? blockops->P_aa(aaHTkey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> a_sector = do_aHT ? sigma->sector(aHTkey) : nullptr;
  shared_ptr<const RASDeterminants> a_det = do_aHT ? a_sector->det() : sigma->space()->det(cc_sector->det()->nelea()-1,cc_sector->det()->neleb());

  const BlockInfo astate(aHTkey.nelea, aHTkey.neleb, nccstates);
  RASBlockVectors sector_r(a_det, astate);

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::AnnihilateAlpha}, {r});
    if (do_aHT) {
      dgemm_("N", "T", a_sector->ndim(), a_sector->mdim(), sector_r.mdim(), 1.0, sector_r.data(), sector_r.ndim(),
                                                &(*S)(0,0,r), S->extent(0), 1.0, a_sector->data(), a_sector->ndim());

      // plus some stuff on Da and Db
    }
#if 0
    if (do_double) {
      // some extra stuff
    }
#endif
  }
}

// should only enter this code if at least one of the terms will be computed
void FormSigmaProdRAS::branch_4(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  // S_beta
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey aHTkey(cckey.nelea+1, cckey.neleb);
  const BlockKey bHTkey(cckey.nelea, cckey.neleb+1);
  const BlockKey bbHTkey(cckey.nelea, cckey.neleb+2);
  const BlockKey abHTkey(cckey.nelea+1, cckey.neleb+1);

  const bool do_aHT = sigma->left()->contains(aHTkey);
  const bool do_bHT = sigma->left()->contains(bHTkey);
  const bool do_bbHT = sigma->left()->contains(bbHTkey);
  const bool do_abHT = sigma->left()->contains(abHTkey);
  assert(do_aHT || do_bHT || do_bbHT || do_abHT);

  shared_ptr<const btas::Tensor3<double>> S = do_bHT ? blockops->S_b(cckey) : nullptr;
  shared_ptr<const btas::TensorN<double,5>> D_b = do_bHT ? blockops->D_b(cckey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_bb = do_bbHT ? blockops->P_bb(bbHTkey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_ab = do_abHT ? blockops->P_ab(abHTkey) : nullptr;
  shared_ptr<const btas::TensorN<double,5>> D_a = do_aHT ? blockops->D_a(cckey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> b_sector = do_bHT ? sigma->sector(bHTkey) : nullptr;
  shared_ptr<const RASDeterminants> b_det = do_bHT ? b_sector->det() : sigma->space()->det(cc_sector->det()->nelea(), cc_sector->det()->neleb()-1);

  const BlockInfo bstate(bHTkey.nelea, bHTkey.neleb, nccstates);
  RASBlockVectors sector_r(b_det, bstate);

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::AnnihilateBeta}, {r});
    if (do_bHT) {
      dgemm_("N", "T", b_sector->ndim(), b_sector->mdim(), sector_r.mdim(), 1.0, sector_r.data(), sector_r.ndim(),
                                                &(*S)(0,0,r), S->extent(0), 1.0, b_sector->data(), b_sector->ndim());
      // plus some displacement stuff
    }
#if 0
    if (do_aHT) {
      // some extra stuff
    }

    if (do_bbHT) {
      // some extra stuff
    }

    if (do_abHT) {
      // some extra stuff
    }
#endif
  }
}


void FormSigmaProdRAS::branch_5(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const BlockKey cckey = cc_sector->left_state().key();
  const int nccstates = cc_sector->mdim();
  RASBlockVectors sector_rs(cc_sector->det(), cc_sector->left_state());
  shared_ptr<const btas::Tensor4<double>> Qaa = blockops->Q_aa(cckey);

  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(cckey);

  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger s)_alpha to the Civecs
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateAlpha,GammaSQ::CreateAlpha}, {s,r});
      dgemm_("N","T", sigma_sector->ndim(), sigma_sector->mdim(), sigma_sector->mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qaa)(0,0,s,r), Qaa->extent(0),
                                                                                        1.0, sigma_sector->data(), sigma_sector->ndim());

      // Then possibly some ET stuff
    }
  }
}

void FormSigmaProdRAS::branch_6(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const BlockKey cckey = cc_sector->left_state().key();
  RASBlockVectors sector_rs(cc_sector->det(), cc_sector->left_state());
  const int nccstates = cc_sector->mdim();

  shared_ptr<const btas::Tensor4<double>> Qbb = blockops->Q_bb(cckey);

  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(cckey);

  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger s)_beta to the Civecs
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateBeta,GammaSQ::CreateBeta}, {s,r});
      dgemm_("N","T", sigma_sector->ndim(), sigma_sector->mdim(), sigma_sector->mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qbb)(0,0,s,r), Qbb->extent(0),
                                                                                        1.0, sigma_sector->data(), sigma_sector->ndim());

      // Then possibly some ET stuff
    }
  }

}

void FormSigmaProdRAS::branch_7(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const int nccstates = cc_sector->mdim();
  const BlockKey cckey = cc_sector->left_state().key();
  const BlockKey flipkey(cckey.nelea-1, cckey.neleb+1);
  const BlockKey bETkey(cckey.nelea-1, cckey.neleb);

  // figure out which sections to do
  const bool do_flip = sigma->left()->contains(flipkey);
  const bool do_bET = sigma->left()->contains(bETkey);
  assert(do_flip || do_bET);

  shared_ptr<const RASDeterminants> flipdet = do_flip ? sigma->sector(flipkey)->det() : sigma->space()->det(cc_sector->det()->nelea()+1,cc_sector->det()->neleb()-1);
  shared_ptr<const btas::Tensor4<double>> Qab = do_flip ? blockops->Q_ab(cckey) : nullptr;

  shared_ptr<RASBlockVectors> sigma_sector = do_flip ? sigma->sector(flipkey) : nullptr;

  RASBlockVectors sector_rs(flipdet, BlockInfo(flipkey.nelea, flipkey.neleb, nccstates));
  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger_alpha s_beta) to the civec
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateBeta,GammaSQ::CreateAlpha}, {s, r});
      if (do_flip)
        dgemm_("N", "T", sigma_sector->ndim(), sigma_sector->mdim(), sector_rs.mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qab)(0,0,s,r), Qab->extent(0),
                                                                                       1.0, sigma_sector->data(), sigma_sector->ndim());
#if 0
      if (do_bET) {
        // some stuff
      }
#endif
    }
  }
}

void FormSigmaProdRAS::branch_8(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const int nccstates = cc_sector->mdim();
  const BlockKey cckey = cc_sector->left_state().key();
  const BlockKey flipkey(cckey.nelea+1, cckey.neleb-1);
  const BlockKey bETkey(cckey.nelea, cckey.neleb-1);

  // figure out which sections to do
  const bool do_flip = sigma->left()->contains(flipkey);
  const bool do_bET = sigma->left()->contains(bETkey);
  assert(do_flip || do_bET);

  shared_ptr<const RASDeterminants> flipdet = do_flip ? sigma->sector(flipkey)->det() : sigma->space()->det(cc_sector->det()->nelea()-1,cc_sector->det()->neleb()+1);
  shared_ptr<const btas::Tensor4<double>> Qab = do_flip ? blockops->Q_ab(flipkey) : nullptr;

  shared_ptr<RASBlockVectors> sigma_sector = do_flip ? sigma->sector(flipkey) : nullptr;

  RASBlockVectors sector_rs(flipdet, BlockInfo(flipkey.nelea, flipkey.neleb, nccstates));
  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger_beta s_alpha) to the civec
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateAlpha,GammaSQ::CreateBeta}, {s, r});
      if (do_flip)
        dgemm_("N", "N", sigma_sector->ndim(), sigma_sector->mdim(), sector_rs.mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qab)(0,0,r,s), Qab->extent(0),
                                                                                       1.0, sigma_sector->data(), sigma_sector->ndim());
#if 0
      if (do_aET) {
        // some stuff
      }
#endif
    }
  }
}

