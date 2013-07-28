//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <utility>

#include <src/meh/gamma_forest.h>

using namespace bagel;
using namespace std;

GammaTree::GammaTree(shared_ptr<const Dvec> ket) : ket_(ket) {
  base_ = make_shared<GammaBranch>();
  const int nops = 4;

  for (int i = 0; i < nops; ++i) {
    base_->branch(i) = make_shared<GammaBranch>();
    for (int j = 0; j < nops; ++j) {
      base_->branch(i)->branch(j) = make_shared<GammaBranch>();
      for (int k = 0; k < nops; ++k) {
        base_->branch(i)->branch(j)->branch(k) = make_shared<GammaBranch>();
      }
    }
  }
}

void GammaTree::compute() {
  const int nops = 4;

  const int nA = ket_->ij();
  const int norb = ket_->det()->norb();

  // Allocation sweep
  for (int i = 0; i < nops; ++i) {
    auto first = base_->branch(i);
    if (!first->active()) continue;
    for (auto& ibra : first->bras()) {
      const int nAp = ibra.second->ij();
      const int nstates = nA * nAp;
      first->gammas().insert(make_pair(ibra.first, make_shared<Matrix>(nstates, norb)));
    }

    for (int j = 0; j < nops; ++j) {
      auto second = first->branch(j);
      if (!second->active()) continue;
      for (auto& jbra : second->bras()) {
        const int nAp = jbra.second->ij();
        const int nstates = nA * nAp;
        second->gammas().insert(make_pair(jbra.first, make_shared<Matrix>(nstates, norb * norb)));
      }

      for (int k = 0; k < nops; ++k) {
        auto third = second->branch(k);
        if (!third->active()) continue;
        for (auto& kbra : third->bras()) {
          const int nAp = kbra.second->ij();
          const int nstates = nA * nAp;
          third->gammas().insert(make_pair(kbra.first, make_shared<Matrix>(nstates, norb * norb * norb)));
        }
      }
    }
  }

  // Computation sweep
  for (int i = 0; i < nops; ++i) {
    auto first = base_->branch(i);
    if (!first->active()) continue;

    for (int a = 0; a < norb; ++a) {
      shared_ptr<const Dvec> avec = apply(ket_, GammaSQ(i), a);
      for (auto& ibra : first->bras()) {
        const int nstates = avec->ij() * ibra.second->ij();
        unique_ptr<double[]> tmp = ddot(ibra.second, avec);
        double* target = first->gammas().find(ibra.first)->second->element_ptr(0,a);
        copy_n(tmp.get(), nstates, target);
      }

      for (int j = 0; j < nops; ++j) {
        auto second = first->branch(j);
        if (!second->active()) continue;

        for (int b = 0; b < norb; ++b) {
          shared_ptr<const Dvec> bvec = apply(avec, GammaSQ(j), b);
          for (auto& jbra : second->bras()) {
            const int nstates = bvec->ij() * jbra.second->ij();
            unique_ptr<double[]> tmp = ddot(jbra.second, bvec);
            double* target = second->gammas().find(jbra.first)->second->element_ptr(0, a + norb*b);
            copy_n(tmp.get(), nstates, target);
          }

          for (int k = 0; k < nops; ++k) {
            auto third = second->branch(k);
            if (!third->active()) continue;

            for (int c = 0; c < norb; ++c) {
              shared_ptr<const Dvec> cvec = apply(bvec, GammaSQ(k), c);
              for (auto& kbra : third->bras()) {
                const int nstates = cvec->ij() * kbra.second->ij();
                unique_ptr<double[]> tmp = ddot(kbra.second, cvec);
                double* target = third->gammas().find(kbra.first)->second->element_ptr(0, a + norb * b + norb * norb * c);
                copy_n(tmp.get(), nstates, target);
              }
            }
          }
        }
      }
    }
  }
}

shared_ptr<const Dvec> GammaTree::apply(shared_ptr<const Dvec> kets, const GammaSQ operation, const int orbital) const {
  const int action = ( (operation==GammaSQ::CreateAlpha || operation==GammaSQ::CreateBeta) ? Create : Annihilate);
  const int spin = ( (operation==GammaSQ::CreateAlpha || operation==GammaSQ::AnnihilateAlpha) ? Alpha : Beta );

  shared_ptr<const Determinants> source_det = kets->det();
  const int norb = source_det->norb();
  const int nstates = kets->ij();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  if (spin == Alpha) {
    shared_ptr<const Determinants> target_det = ( action == Annihilate ? source_det->remalpha() : source_det->addalpha() );

    auto out = make_shared<Dvec>(target_det, nstates);

    const int target_lena = target_det->lena();
    const int target_lenb = target_det->lenb();

    for (int i = 0; i < nstates; ++i) {
      double* target_base = out->data(i)->data();
      const double* source_base = kets->data(i)->data();

      for (auto& iter : (action == Annihilate ? source_det->phidowna(orbital) : source_det->phiupa(orbital)) ) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + target_lenb * iter.target;
        const double* source = source_base + source_lenb * iter.source;
        daxpy_(target_lenb, sign, source, 1, target, 1);
      }
    }

    return out;
  }
  else {
    shared_ptr<const Determinants> target_det = (action == Annihilate ? source_det->rembeta() : source_det->addbeta());

    const int target_lena = target_det->lena();
    const int target_lenb = target_det->lenb();

    auto out = make_shared<Dvec>(target_det, nstates);

    for (int i = 0; i < target_lena; ++i) {
      for (int istate = 0; istate < nstates; ++istate) {
        double* target_base = out->data(istate)->element_ptr(0,i);
        const double* source_base = kets->data(istate)->element_ptr(0,i);
        for (auto& iter : (action == Annihilate ? source_det->phidownb(orbital) : source_det->phiupb(orbital))) {
          const double sign = static_cast<double>(iter.sign);
          target_base[iter.target] += sign * source_base[iter.source];
        }
      }
    }

    return out;
  }
}

unique_ptr<double[]> GammaTree::ddot(shared_ptr<const Dvec> bras, shared_ptr<const Dvec> kets) const {
  const int nbras = bras->ij();
  const int nkets = kets->ij();

  unique_ptr<double[]> out(new double[nbras*nkets]);
  double* odata = out.get();

  for (int iket = 0; iket < nkets; ++iket) {
    for (int jbra = 0; jbra < nbras; ++jbra, ++odata) {
      *odata = bras->data(jbra)->ddot(*kets->data(iket));
    }
  }

  return out;
}
