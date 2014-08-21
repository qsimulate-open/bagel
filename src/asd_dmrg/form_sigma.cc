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

#include <src/math/sparsematrix.h>
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
      aET_branch(cc_sector, sigma, blockops);
      ptime.tick_print("aET-branch");
    }

    if (do_bET || do_bbET || do_abET) {
      bET_branch(cc_sector, sigma, blockops);
      ptime.tick_print("bET-branch");
    }

    if (do_aHT || do_aaHT) {
      aHT_branch(cc_sector, sigma, blockops);
      ptime.tick_print("aHT-branch");
    }

    if (do_bHT || do_bbHT || do_abHT) {
      bHT_branch(cc_sector, sigma, blockops);
      ptime.tick_print("bHT-branch");
    }

    // always compute these
    aexc_branch(cc_sector, sigma, blockops);
    bexc_branch(cc_sector, sigma, blockops);
    ptime.tick_print("exc-branches");

    if (do_flipup) {
      abflip_branch(cc_sector, sigma, blockops);
      ptime.tick_print("abflip-branch");
    }

    if (do_flipdn) {
      baflip_branch(cc_sector, sigma, blockops);
      ptime.tick_print("baflip-branch");
    }

    if (do_aET) {
      compute_sigma_3aET(cc_sector, sigma, jop);
      ptime.tick_print("sigma-3aET");
    }

    if (do_aHT) {
      compute_sigma_3aHT(cc_sector, sigma, jop);
      ptime.tick_print("sigma-3aHT");
    }

    if (do_bET) {
      compute_sigma_3bET(cc_sector, sigma, jop);
      ptime.tick_print("sigma-3bET");
    }

    if (do_bHT) {
      compute_sigma_3bHT(cc_sector, sigma, jop);
      ptime.tick_print("sigma-3bHT");
    }
  }
}

void FormSigmaProdRAS::aET_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
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

  shared_ptr<RASBlockVectors> double_sector = do_double ? sigma->sector(doubleETkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_double = do_double ? make_shared<RASBlockVectors>(double_sector->det(), BlockInfo(doubleETkey.nelea, doubleETkey.neleb, nccstates))
                                                     : nullptr;

  const BlockInfo singlestate(singleETkey.nelea, singleETkey.neleb, nccstates);
  RASBlockVectors sector_r(single_det, BlockInfo(singleETkey.nelea, singleETkey.neleb, nccstates));

  const int phase = (1 - (((cc_sector->det()->nelea()+cc_sector->det()->neleb())%2) << 1));

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::CreateAlpha}, {r});
    if (do_single)
      dgemm_("N", "N", single_sector->ndim(), single_sector->mdim(), sector_r.mdim(), phase, sector_r.data(), sector_r.ndim(),
                                                            &(*S)(0,0,r), S->extent(0), 1.0, single_sector->data(), single_sector->ndim());
    if (do_double) {
      for (int s = 0; s < r; ++s) {
        tmp_double->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_double->civec(ist), {GammaSQ::CreateAlpha}, {s});
        Matrix pa(P->extent(0), P->extent(1));
        blas::ax_plus_y_n(1.0, &(*P)(0,0,r,s), pa.size(), pa.data());
        blas::ax_plus_y_n(-1.0, &(*P)(0,0,s,r), pa.size(), pa.data());
        dgemm_("N", "T", double_sector->ndim(), double_sector->mdim(), nccstates, 1.0, tmp_double->data(), tmp_double->ndim(),
                                                            pa.data(), pa.ndim(), 1.0, double_sector->data(), double_sector->ndim());
      }
    }
  }
}

void FormSigmaProdRAS::bET_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
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
  shared_ptr<const RASDeterminants> b_det = do_b ? b_sector->det() : sigma->space()->det(cc_sector->det()->nelea(), cc_sector->det()->neleb()+1);

  shared_ptr<RASBlockVectors> bb_sector = do_bb ? sigma->sector(bbETkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_bb = do_bb ? make_shared<RASBlockVectors>(bb_sector->det(), BlockInfo(bbETkey.nelea, bbETkey.neleb, nccstates)) : nullptr;

  shared_ptr<RASBlockVectors> ab_sector = do_ab ? sigma->sector(abETkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_ab = do_ab ? make_shared<RASBlockVectors>(ab_sector->det(), BlockInfo(abETkey.nelea, abETkey.neleb, nccstates)) : nullptr;

  const BlockInfo bstate(bETkey.nelea, bETkey.neleb, nccstates);
  RASBlockVectors sector_r(b_det, BlockInfo(bETkey.nelea, bETkey.neleb, nccstates));

  const int phase = (1 - (((cc_sector->det()->nelea()+cc_sector->det()->neleb())%2) << 1));

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::CreateBeta}, {r});
    if (do_b)
      dgemm_("N", "N", b_sector->ndim(), b_sector->mdim(), sector_r.mdim(), phase, sector_r.data(), sector_r.ndim(),
                                                  &(*S)(0,0,r), S->extent(0), 1.0, b_sector->data(), b_sector->ndim());
    if (do_bb) {
      for (int s = 0; s < r; ++s) {
        tmp_bb->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_bb->civec(ist), {GammaSQ::CreateBeta}, {s});
        Matrix pbb(P_bb->extent(0), P_bb->extent(1));
        blas::ax_plus_y_n(1.0, &(*P_bb)(0,0,r,s), pbb.size(), pbb.data());
        blas::ax_plus_y_n(-1.0, &(*P_bb)(0,0,s,r), pbb.size(), pbb.data());
        dgemm_("N", "T", bb_sector->ndim(), bb_sector->mdim(), nccstates, 1.0, tmp_bb->data(), tmp_bb->ndim(),
                                                            pbb.data(), pbb.ndim(), 1.0, bb_sector->data(), bb_sector->ndim());
      }
    }

    if (do_ab) {
      for (int s = 0; s < rnorb; ++s) {
        tmp_ab->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_ab->civec(ist), {GammaSQ::CreateAlpha}, {s});
        dgemm_("N", "T", ab_sector->ndim(), ab_sector->mdim(), nccstates, 1.0, tmp_ab->data(), tmp_ab->ndim(),
                                      &(*P_ab)(0,0,r,s), P_ab->extent(0), 1.0, ab_sector->data(), ab_sector->ndim());
      }
    }
  }
}

void FormSigmaProdRAS::aHT_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  // S_alpha
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey aHTkey(cckey.nelea+1, cckey.neleb);
  const BlockKey aaHTkey(cckey.nelea+2, cckey.neleb);

  const bool do_aHT = sigma->left()->contains(aHTkey);
  const bool do_aaHT = sigma->left()->contains(aaHTkey);
  // not sure how you could do a aa but not a a, but that's a different problem
  assert(do_aHT || do_aaHT);

  shared_ptr<const btas::Tensor3<double>> S = do_aHT ? blockops->S_a(cckey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P = do_aaHT ? blockops->P_aa(aaHTkey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> a_sector = do_aHT ? sigma->sector(aHTkey) : nullptr;
  shared_ptr<const RASDeterminants> a_det = do_aHT ? a_sector->det() : sigma->space()->det(cc_sector->det()->nelea()-1,cc_sector->det()->neleb());

  shared_ptr<RASBlockVectors> aa_sector = do_aaHT ? sigma->sector(aaHTkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_aa = do_aaHT ? make_shared<RASBlockVectors>(aa_sector->det(), BlockInfo(aaHTkey.nelea, aaHTkey.neleb, nccstates)) : nullptr;

  RASBlockVectors sector_r(a_det, BlockInfo(aHTkey.nelea, aHTkey.neleb, nccstates));

  const int phase = (1 - (((sector_r.det()->nelea()+sector_r.det()->neleb())%2) << 1));

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::AnnihilateAlpha}, {r});
    if (do_aHT) {
      dgemm_("N", "T", a_sector->ndim(), a_sector->mdim(), sector_r.mdim(), phase, sector_r.data(), sector_r.ndim(),
                                                  &(*S)(0,0,r), S->extent(0), 1.0, a_sector->data(), a_sector->ndim());
    }

    if (do_aaHT) {
      for (int s = 0; s < r; ++s) {
        tmp_aa->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_aa->civec(ist), {GammaSQ::AnnihilateAlpha}, {s});
        Matrix pa(P->extent(0), P->extent(1));
        blas::ax_plus_y_n(1.0, &(*P)(0,0,s,r), pa.size(), pa.data());
        blas::ax_plus_y_n(-1.0, &(*P)(0,0,r,s), pa.size(), pa.data());
        dgemm_("N", "N", aa_sector->ndim(), aa_sector->mdim(), nccstates, 1.0, tmp_aa->data(), tmp_aa->ndim(),
                                                    pa.data(), pa.ndim(), 1.0, aa_sector->data(), aa_sector->ndim());
      }
    }
  }
}

// should only enter this code if at least one of the terms will be computed
void FormSigmaProdRAS::bHT_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  // S_beta
  const BlockKey cckey = cc_sector->left_state().key();

  const BlockKey bHTkey(cckey.nelea, cckey.neleb+1);
  const BlockKey bbHTkey(cckey.nelea, cckey.neleb+2);
  const BlockKey abHTkey(cckey.nelea+1, cckey.neleb+1);

  const bool do_bHT = sigma->left()->contains(bHTkey);
  const bool do_bbHT = sigma->left()->contains(bbHTkey);
  const bool do_abHT = sigma->left()->contains(abHTkey);
  assert(do_bHT || do_bbHT || do_abHT);

  shared_ptr<const btas::Tensor3<double>> S = do_bHT ? blockops->S_b(cckey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_bb = do_bbHT ? blockops->P_bb(bbHTkey) : nullptr;
  shared_ptr<const btas::Tensor4<double>> P_ab = do_abHT ? blockops->P_ab(abHTkey) : nullptr;

  const int nccstates = cc_sector->nstates();

  shared_ptr<RASBlockVectors> b_sector = do_bHT ? sigma->sector(bHTkey) : nullptr;
  shared_ptr<const RASDeterminants> b_det = do_bHT ? b_sector->det() : sigma->space()->det(cc_sector->det()->nelea(), cc_sector->det()->neleb()-1);

  RASBlockVectors sector_r(b_det, BlockInfo(bHTkey.nelea, bHTkey.neleb, nccstates));

  shared_ptr<RASBlockVectors> bb_sector = do_bbHT ? sigma->sector(bbHTkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_bb = do_bbHT ? make_shared<RASBlockVectors>(bb_sector->det(), BlockInfo(bbHTkey.nelea, bbHTkey.neleb, nccstates)) : nullptr;

  shared_ptr<RASBlockVectors> ab_sector = do_abHT ? sigma->sector(abHTkey) : nullptr;
  shared_ptr<RASBlockVectors> tmp_ab = do_abHT ? make_shared<RASBlockVectors>(ab_sector->det(), BlockInfo(abHTkey.nelea, abHTkey.neleb, nccstates)) : nullptr;

  const int phase = (1 - (((sector_r.det()->nelea()+sector_r.det()->neleb())%2) << 1));

  for (int r = 0; r < rnorb; ++r) {
    sector_r.zero();
    for (int ist = 0; ist < nccstates; ++ist)
      apply(1.0, cc_sector->civec(ist), sector_r.civec(ist), {GammaSQ::AnnihilateBeta}, {r});
    if (do_bHT) {
      dgemm_("N", "T", b_sector->ndim(), b_sector->mdim(), sector_r.mdim(), phase, sector_r.data(), sector_r.ndim(),
                                                  &(*S)(0,0,r), S->extent(0), 1.0, b_sector->data(), b_sector->ndim());
    }

    if (do_bbHT) {
      for (int s = 0; s < r; ++s) {
        tmp_bb->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_bb->civec(ist), {GammaSQ::AnnihilateBeta}, {s});
        Matrix pbb(P_bb->extent(0), P_bb->extent(1));
        blas::ax_plus_y_n(1.0, &(*P_bb)(0,0,s,r), pbb.size(), pbb.data());
        blas::ax_plus_y_n(-1.0, &(*P_bb)(0,0,r,s), pbb.size(), pbb.data());
        dgemm_("N", "N", bb_sector->ndim(), bb_sector->mdim(), nccstates, 1.0, tmp_bb->data(), tmp_bb->ndim(),
                                                  pbb.data(), pbb.ndim(), 1.0, bb_sector->data(), bb_sector->ndim());
      }
    }

    if (do_abHT) {
      for (int s = 0; s < rnorb; ++s) {
        tmp_ab->zero();
        for (int ist = 0; ist < nccstates; ++ist)
          apply(1.0, sector_r.civec(ist), tmp_ab->civec(ist), {GammaSQ::AnnihilateAlpha}, {s});
        dgemm_("N", "N", ab_sector->ndim(), ab_sector->mdim(), nccstates, -1.0, tmp_ab->data(), tmp_ab->ndim(),
                                       &(*P_ab)(0,0,r,s), P_ab->extent(0), 1.0, ab_sector->data(), ab_sector->ndim());
      }
    }
  }
}


void FormSigmaProdRAS::aexc_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
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

    }
  }
}

void FormSigmaProdRAS::bexc_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
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

    }
  }

}

void FormSigmaProdRAS::abflip_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const int nccstates = cc_sector->mdim();
  const BlockKey cckey = cc_sector->left_state().key();
  const BlockKey flipkey(cckey.nelea-1, cckey.neleb+1);

  assert(sigma->left()->contains(flipkey));

  shared_ptr<const RASDeterminants> flipdet = sigma->sector(flipkey)->det();
  shared_ptr<const btas::Tensor4<double>> Qab = blockops->Q_ab(cckey);

  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(flipkey);

  RASBlockVectors sector_rs(flipdet, BlockInfo(flipkey.nelea, flipkey.neleb, nccstates));
  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger_alpha s_beta) to the civec
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateBeta,GammaSQ::CreateAlpha}, {s, r});
      dgemm_("N", "T", sigma_sector->ndim(), sigma_sector->mdim(), sector_rs.mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qab)(0,0,s,r), Qab->extent(0),
                                                                                     1.0, sigma_sector->data(), sigma_sector->ndim());
    }
  }
}

void FormSigmaProdRAS::baflip_branch(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const BlockOperators> blockops) const {
  ApplyOperator apply;
  const int rnorb = cc_sector->det()->norb();

  const int nccstates = cc_sector->mdim();
  const BlockKey cckey = cc_sector->left_state().key();
  const BlockKey flipkey(cckey.nelea+1, cckey.neleb-1);
  assert(sigma->left()->contains(flipkey));

  shared_ptr<const RASDeterminants> flipdet = sigma->sector(flipkey)->det();
  shared_ptr<const btas::Tensor4<double>> Qab = blockops->Q_ab(flipkey);

  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(flipkey);

  RASBlockVectors sector_rs(flipdet, BlockInfo(flipkey.nelea, flipkey.neleb, nccstates));
  for (int r = 0; r < rnorb; ++r) {
    for (int s = 0; s < rnorb; ++s) {
      // apply (r^dagger_beta s_alpha) to the civec
      sector_rs.zero();
      for (int ist = 0; ist < nccstates; ++ist)
        apply(1.0, cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateAlpha,GammaSQ::CreateBeta}, {s, r});

      dgemm_("N", "N", sigma_sector->ndim(), sigma_sector->mdim(), sector_rs.mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qab)(0,0,r,s), Qab->extent(0),
                                                                                     1.0, sigma_sector->data(), sigma_sector->ndim());
    }
  }
}

void FormSigmaProdRAS::compute_sigma_3aET(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const DimerJop> jop) const {
  const BlockKey aETkey(cc_sector->left_state().nelea-1, cc_sector->left_state().neleb);
  assert(sigma->left()->contains(aETkey));
  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(aETkey);
  const int nccstates = cc_sector->mdim();
  const BlockInfo tmpinfo(aETkey.nelea, aETkey.neleb, nccstates);
  RASBlockVectors tmp_sector(sigma_sector->det(), tmpinfo);

  const int lnorb = jop->monomer_jop<1>()->nocc();
  const int rnorb = jop->monomer_jop<0>()->nocc();

  shared_ptr<const Matrix> J = jop->coulomb_matrix<0,1,0,0>();
  shared_ptr<const btas::Tensor3<double>> gamma_a = sigma->left()->coupling({GammaSQ::CreateAlpha}).at({cc_sector->left_state(), aETkey}).data;

  assert(gamma_a->extent(0)==nccstates && gamma_a->extent(1)==sigma_sector->mdim());

  const int phase = (1 - (((cc_sector->det()->nelea()+cc_sector->det()->neleb())%2) << 1));

  for (int p = 0; p < lnorb; ++p ) {
    tmp_sector.zero();
    auto Jp = make_shared<btas::Tensor3<double>>(rnorb, rnorb, rnorb);
    copy_n(J->element_ptr(0, p), rnorb*rnorb*rnorb, Jp->data());
    for (int ist = 0; ist < nccstates; ++ist) {
      resolve_S_adag_adag_a(cc_sector->civec(ist), tmp_sector.civec(ist), Jp);
      resolve_S_abb(cc_sector->civec(ist), tmp_sector.civec(ist), Jp);
    }

    dgemm_("N", "N", sigma_sector->ndim(), sigma_sector->mdim(), cc_sector->mdim(), phase, tmp_sector.data(), tmp_sector.ndim(),
                                              &(*gamma_a)(0,0,p), gamma_a->extent(0), 1.0, sigma_sector->data(), sigma_sector->ndim());
  }
}

void FormSigmaProdRAS::compute_sigma_3aHT(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const DimerJop> jop) const {
  const BlockKey aHTkey(cc_sector->left_state().nelea+1, cc_sector->left_state().neleb);
  assert(sigma->left()->contains(aHTkey));
  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(aHTkey);
  const int nccstates = cc_sector->mdim();
  const BlockInfo tmpinfo(aHTkey.nelea, aHTkey.neleb, nccstates);
  RASBlockVectors tmp_sector(sigma_sector->det(), tmpinfo);

  const int lnorb = jop->monomer_jop<1>()->nocc();
  const int rnorb = jop->monomer_jop<0>()->nocc();

  shared_ptr<const Matrix> J = jop->coulomb_matrix<0,1,0,0>();
  shared_ptr<const btas::Tensor3<double>> gamma_a = sigma->left()->coupling({GammaSQ::CreateAlpha}).at({aHTkey, cc_sector->left_state()}).data;

  assert(gamma_a->extent(1)==nccstates && gamma_a->extent(0)==sigma_sector->mdim());

  const int phase = (1 - (((tmp_sector.det()->nelea()+tmp_sector.det()->neleb())%2) << 1));

  for (int p = 0; p < lnorb; ++p ) {
    tmp_sector.zero();
    auto Jp = make_shared<btas::Tensor3<double>>(rnorb, rnorb, rnorb);
    copy_n(J->element_ptr(0, p), rnorb*rnorb*rnorb, Jp->data());
    for (int ist = 0; ist < nccstates; ++ist) {
      resolve_S_adag_a_a(cc_sector->civec(ist), tmp_sector.civec(ist), Jp);
      resolve_S_abb(cc_sector->civec(ist), tmp_sector.civec(ist), Jp);
    }

    dgemm_("N", "T", sigma_sector->ndim(), sigma_sector->mdim(), cc_sector->mdim(), phase, tmp_sector.data(), tmp_sector.ndim(),
                                              &(*gamma_a)(0,0,p), gamma_a->extent(0), 1.0, sigma_sector->data(), sigma_sector->ndim());
  }
}

void FormSigmaProdRAS::compute_sigma_3bET(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const DimerJop> jop) const {
  const BlockKey bETkey(cc_sector->left_state().nelea, cc_sector->left_state().neleb-1);
  assert(sigma->left()->contains(bETkey));
  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(bETkey);
  const int nccstates = cc_sector->mdim();
  const BlockInfo tmpinfo(bETkey.nelea, bETkey.neleb, nccstates);
  RASBlockVectors tmp_sector(sigma_sector->det(), tmpinfo);

  const int lnorb = jop->monomer_jop<1>()->nocc();
  const int rnorb = jop->monomer_jop<0>()->nocc();

  shared_ptr<const Matrix> J = jop->coulomb_matrix<0,1,0,0>();
  shared_ptr<const btas::Tensor3<double>> gamma_b = sigma->left()->coupling({GammaSQ::CreateBeta}).at({cc_sector->left_state(), bETkey}).data;

  const int phase = (1 - (((cc_sector->det()->nelea()+cc_sector->det()->neleb())%2) << 1));

  shared_ptr<RASCivec> tmp_trans = tmp_sector.civec(0).transpose();

  for (int p = 0; p < lnorb; ++p ) {
    tmp_sector.zero();
    auto Jp = make_shared<btas::Tensor3<double>>(rnorb, rnorb, rnorb);
    copy_n(J->element_ptr(0, p), rnorb*rnorb*rnorb, Jp->data());
    for (int ist = 0; ist < nccstates; ++ist) {
      auto cc_trans = cc_sector->civec(ist).transpose();
      tmp_trans->zero();
      resolve_S_adag_adag_a(*cc_trans, *tmp_trans, Jp);
      resolve_S_abb(*cc_trans, *tmp_trans, Jp);
      tmp_sector.civec(ist).ax_plus_y(1.0, *tmp_trans->transpose());
    }

    dgemm_("N", "N", sigma_sector->ndim(), sigma_sector->mdim(), cc_sector->mdim(), phase, tmp_sector.data(), tmp_sector.ndim(),
                                              &(*gamma_b)(0,0,p), gamma_b->extent(0), 1.0, sigma_sector->data(), sigma_sector->ndim());
  }
}


void FormSigmaProdRAS::compute_sigma_3bHT(shared_ptr<const RASBlockVectors> cc_sector, shared_ptr<ProductRASCivec> sigma, shared_ptr<const DimerJop> jop) const {
  const BlockKey bHTkey(cc_sector->left_state().nelea, cc_sector->left_state().neleb+1);
  assert(sigma->left()->contains(bHTkey));
  shared_ptr<RASBlockVectors> sigma_sector = sigma->sector(bHTkey);
  const int nccstates = cc_sector->mdim();
  const BlockInfo tmpinfo(bHTkey.nelea, bHTkey.neleb, nccstates);
  RASBlockVectors tmp_sector(sigma_sector->det(), tmpinfo);

  const int lnorb = jop->monomer_jop<1>()->nocc();
  const int rnorb = jop->monomer_jop<0>()->nocc();

  shared_ptr<const Matrix> J = jop->coulomb_matrix<0,1,0,0>();
  shared_ptr<const btas::Tensor3<double>> gamma_b = sigma->left()->coupling({GammaSQ::CreateBeta}).at({bHTkey, cc_sector->left_state()}).data;

  const int phase = (1 - (((tmp_sector.det()->nelea()+tmp_sector.det()->neleb())%2) << 1));

  shared_ptr<RASCivec> tmp_trans = tmp_sector.civec(0).transpose();

  for (int p = 0; p < lnorb; ++p ) {
    tmp_sector.zero();
    auto Jp = make_shared<btas::Tensor3<double>>(rnorb, rnorb, rnorb);
    copy_n(J->element_ptr(0, p), rnorb*rnorb*rnorb, Jp->data());
    for (int ist = 0; ist < nccstates; ++ist) {
      auto cc_trans = cc_sector->civec(ist).transpose();
      tmp_trans->zero();
      resolve_S_adag_a_a(*cc_trans, *tmp_trans, Jp);
      resolve_S_abb(*cc_trans, *tmp_trans, Jp);
      tmp_sector.civec(ist).ax_plus_y(1.0, *tmp_trans->transpose());
    }

    dgemm_("N", "T", sigma_sector->ndim(), sigma_sector->mdim(), cc_sector->mdim(), phase, tmp_sector.data(), tmp_sector.ndim(),
                                              &(*gamma_b)(0,0,p), gamma_b->extent(0), 1.0, sigma_sector->data(), sigma_sector->ndim());
  }
}


// resolves (S_p)_alpha = \sum_{ijk} (pk|ij) (k^+)_alpha i^+_alpha j_alpha
void FormSigmaProdRAS::resolve_S_adag_adag_a(const RASCivecView cc, RASCivecView sigma, shared_ptr<btas::Tensor3<double>> Jp) const {
  shared_ptr<const RASDeterminants> sdet = cc.det();
  shared_ptr<const RASDeterminants> tdet = sigma.det();

  const int norb = sdet->norb();
  assert(norb == tdet->norb());

  const size_t sla = sdet->lena();

  // k^+_alpha i^+_alpha j_alpha portion. the way easy part
  Matrix F(sla, batchsize_);
  for (auto& targetspace : *tdet->stringspacea()) {
    const int nbatches = (targetspace->size()-1)/batchsize_ + 1;
    for (int batch = 0; batch < nbatches; ++batch) {
      const size_t batchstart = batch*batchsize_;
      const size_t batchlength = min(static_cast<size_t>(batchsize_), targetspace->size() - batchstart);

      F.zero();

      // first build an F
      for (size_t ia = 0; ia < batchlength; ++ia) {
        double* const fdata = F.element_ptr(0, ia);
        const bitset<nbit__> tabit = targetspace->strings(ia);
        for (int k = 0; k < norb; ++k) {
          if (!tabit[k]) continue;
          bitset<nbit__> tmpbit = tabit ^ bitset<nbit__>(1 << k);
          if (!sdet->allowed(tmpbit)) continue;
          const int kphase = sign(tmpbit, k);
          const size_t tmpia = sdet->lexical_offset<0>(tmpbit);
          for (auto& iterij : sdet->phia(tmpia)) {
            const int i = iterij.ij%norb;
            const int j = iterij.ij/norb;
            if (i < k)
              fdata[iterij.source] += static_cast<double>(kphase*iterij.sign) * ((*Jp)(k, i, j) - (*Jp)(i, k, j));
          }
        }
      }

      // Sp(beta, alpha) += \sum_alpha' C(beta, alpha') * F(alpha', alpha)
      for (auto& ccblock : cc.blocks()) {
        if (!ccblock) continue;
        if (!tdet->allowed(targetspace, ccblock->stringsb())) continue;
        shared_ptr<RASBlock<double>> target_block = sigma.block(ccblock->stringsb(), targetspace);

        assert(ccblock->lenb() == target_block->lenb());
        dgemm_("N", "N", target_block->lenb(), batchlength, ccblock->lena(), 1.0, ccblock->data(), ccblock->lenb(),
                  F.element_ptr(ccblock->stringsa()->offset(), 0), F.ndim(), 1.0, target_block->data() + batchstart * target_block->lenb(), target_block->lenb());
      }
    }
  }
}

// resolves (S_p)_alpha^\dagger = \sum_{ijk} (pk|ij) (i^+)_alpha j_alpha k_alpha
void FormSigmaProdRAS::resolve_S_adag_a_a(const RASCivecView cc, RASCivecView sigma, shared_ptr<btas::Tensor3<double>> Jp) const {
  shared_ptr<const RASDeterminants> sdet = cc.det();
  shared_ptr<const RASDeterminants> tdet = sigma.det();

  const int norb = sdet->norb();
  assert(norb == tdet->norb());

  const size_t sla = sdet->lena();

  // i^+_alpha j^_alpha k_alpha portion. the way easy part
  Matrix F(sla, batchsize_);
  for (auto& targetspace : *tdet->stringspacea()) {
    const int nbatches = (targetspace->size()-1)/batchsize_ + 1;
    for (int batch = 0; batch < nbatches; ++batch) {
      const size_t batchstart = batch*batchsize_;
      const size_t batchlength = min(static_cast<size_t>(batchsize_), targetspace->size() - batchstart);

      F.zero();

      // first build an F
      const size_t offset = batchstart + targetspace->offset();
      for (size_t ia = 0; ia < batchlength; ++ia) {
        double* const fdata = F.element_ptr(0, ia);
        for (auto& iterij : tdet->phia(ia + offset)) {
          const bitset<nbit__> tmpbit = tdet->string_bits_a(iterij.source);
          const int i = iterij.ij%norb;
          const int j = iterij.ij/norb;
          for (int k = j+1; k < norb; ++k) {
            if (tmpbit[k]) continue;
            const bitset<nbit__> sabit = tmpbit ^ bitset<nbit__>(1 << k);
            if (!sdet->allowed(sabit)) continue;
            const int kphase = sign(sabit, k);
            const size_t source_ia = sdet->lexical_offset<0>(sabit);
            fdata[source_ia] += static_cast<double>(kphase*iterij.sign) * ((*Jp)(k, i, j) - (*Jp)(j, i, k));
          }
        }
      }

      // Sp(beta, alpha) += \sum_alpha' C(beta, alpha') * F(alpha', alpha)
      for (auto& ccblock : cc.blocks()) {
        if (!ccblock) continue;
        if (!tdet->allowed(targetspace, ccblock->stringsb())) continue;
        shared_ptr<RASBlock<double>> target_block = sigma.block(ccblock->stringsb(), targetspace);

        assert(ccblock->lenb() == target_block->lenb());
        dgemm_("N", "N", target_block->lenb(), batchlength, ccblock->lena(), 1.0, ccblock->data(), ccblock->lenb(),
                  F.element_ptr(ccblock->stringsa()->offset(), 0), F.ndim(), 1.0, target_block->data() + batchstart * target_block->lenb(), target_block->lenb());
      }
    }
  }
}

// computes \sum_{ijk} (k)_alpha i^+_beta j_beta (pk|ij)
// (k) can be creation or annihilation and is determined based on the electron count of sigma
void FormSigmaProdRAS::resolve_S_abb(const RASCivecView cc, RASCivecView sigma, shared_ptr<btas::Tensor3<double>> Jp) const {
  shared_ptr<const RASDeterminants> sdet = cc.det();
  shared_ptr<const RASDeterminants> tdet = sigma.det();

  assert(abs(sdet->nelea()-tdet->nelea())==1);
  assert(sdet->neleb()==tdet->neleb());

  const size_t na_target = tdet->nelea();

  const int norb = sdet->norb();
  assert(norb == tdet->norb());

  // k^?_alpha i^+_beta j_beta portion. the harder part
  for (int k = 0; k < norb; ++k) {
    // map<taga, map<tagb, pair<Cprime, list of target destinations>>>
    map<int, map<int, shared_ptr<Matrix>>> Cp_map;
    map<int, vector<tuple<size_t, bitset<nbit__>>>> RI_map;

    // gathering operations
    for (auto& source_aspace : *sdet->stringspacea()) {
      vector<tuple<size_t, int>> kmap;
      for (size_t ia = 0; ia < source_aspace->size(); ++ia) {
        const bitset<nbit__> sabit = sdet->string_bits_a(ia+source_aspace->offset());
        const bitset<nbit__> tabit = sabit ^ bitset<nbit__>(1 << k); // flips occupation k
        // counting nelea should tell me whether flipping k was in the right direction or not
        if (tabit.count()==na_target) {
          // if too many holes or particles
          if (!tdet->allowed(tabit)) continue;

          const int phase = sign(tabit, k);
          const size_t target_lex = tdet->lexical_zero<0>(tabit);
          kmap.emplace_back(ia, phase);
          RI_map[source_aspace->tag()].emplace_back(target_lex, tabit);
        }
      }

      for (auto& iblock : cc.allowed_blocks<0>(source_aspace)) {
        auto cp_block = make_shared<Matrix>(iblock->lenb(), kmap.size());
        for (size_t ia = 0; ia < kmap.size(); ++ia)
          blas::ax_plus_y_n(get<1>(kmap[ia]), iblock->data() + iblock->lenb()*get<0>(kmap[ia]), iblock->lenb(), cp_block->element_ptr(0, ia));
        Cp_map[source_aspace->tag()].emplace(iblock->stringsb()->tag(), cp_block);
      }
    }

    for (auto& target_bspace : *tdet->stringspaceb()) {
      for (auto& aspace_RI : RI_map) {
        vector<int> reduced_positions;
        vector<tuple<size_t,bitset<nbit__>>> reduced_RI;
        int current = 0;
        for (auto& i : aspace_RI.second) {
          shared_ptr<const RASString> target_aspace = tdet->space<0>(get<1>(i));
          if (tdet->allowed(target_bspace, target_aspace)) {
            reduced_positions.push_back(current);
            reduced_RI.push_back(i);
          }
          ++current;
        }

        if (reduced_RI.empty()) continue;

        for (auto& source_bspace : *sdet->stringspaceb()) {
          shared_ptr<const Matrix> full_cp = Cp_map[aspace_RI.first][source_bspace->tag()];
          if (!full_cp) continue;
          auto reduced_cp = make_shared<Matrix>(full_cp->ndim(), reduced_RI.size());
          current = 0;
          for (auto& i : reduced_positions)
            copy_n(full_cp->element_ptr(0,i), full_cp->ndim(), reduced_cp->element_ptr(0,current++));

          // Now build an F matrix in sparse format
          map<pair<int, int>, double> sparse_coords;
          for (size_t ib = 0; ib < target_bspace->size(); ++ib) {
            for (auto& iterij : tdet->phib(ib+target_bspace->offset())) {
              if (iterij.source >= source_bspace->offset() && iterij.source < source_bspace->offset()+source_bspace->size()) {
                const int i = iterij.ij%norb;
                const int j = iterij.ij/norb;
                sparse_coords[{ib,iterij.source-source_bspace->offset()}] += iterij.sign * (*Jp)(k,i,j);
              }
            }
          }

          if (!sparse_coords.empty()) {
            SparseMatrix sparseF(target_bspace->size(), source_bspace->size(), sparse_coords);
            Matrix V(sparseF * *reduced_cp);

            // scatter
            current = 0;
            for (auto& ireduced : reduced_RI) {
              shared_ptr<const RASString> target_aspace = tdet->space<0>(get<1>(ireduced));
              shared_ptr<RASBlock<double>> target_block = sigma.block(target_bspace, target_aspace);
              blas::ax_plus_y_n(1.0, V.element_ptr(0,current++), V.ndim(), target_block->data() + target_block->lenb()*get<0>(ireduced));
            }
          }
        }
      }
    }
  }
}

