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

/* Implementing the method as described by Olsen */
// This is how this is accessed from RASCI
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

    // and now for the interaction part
    shared_ptr<const btas::Tensor4<double>> Qaa = blockops->Q_aa(sector.first);
    shared_ptr<const btas::Tensor4<double>> Qbb = blockops->Q_bb(sector.first);

    RASBlockVectors sector_rs(cc_sector->det(), cc_sector->left_state());
    ApplyOperator apply;

    const int rnorb = jop->monomer_jop<0>()->nocc();

    for (int r = 0; r < rnorb; ++r) {
      for (int s = 0; s < rnorb; ++s) {
        // apply (r^dagger s)_alpha to the Civecs
        sector_rs.zero();
        for (int ist = 0; ist < nsecstates; ++ist)
          apply(cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateAlpha,GammaSQ::CreateAlpha}, {r,s});
        dgemm_("N","T", sigma_sector->ndim(), sigma_sector->mdim(), sigma_sector->mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qaa)(0,0,r,s), nsecstates,
                                                                                          1.0, sigma_sector->data(), sigma_sector->ndim());
        // apply (r^dagger s)_beta to the Civecs
        sector_rs.zero();
        for (int ist = 0; ist < nsecstates; ++ist)
          apply(cc_sector->civec(ist), sector_rs.civec(ist), {GammaSQ::AnnihilateBeta,GammaSQ::CreateBeta}, {r,s});
        dgemm_("N","T", sigma_sector->ndim(), sigma_sector->mdim(), sigma_sector->mdim(), 1.0, sector_rs.data(), sector_rs.ndim(), &(*Qbb)(0,0,r,s), nsecstates,
                                                                                          1.0, sigma_sector->data(), sigma_sector->ndim());
      }
    }
    ptime.tick_print("mixed form sigma");
  }
}
