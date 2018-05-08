//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relreference.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/mat1e/mixedbasis.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/mat1e/giao/reloverlap_london.h>
#include <src/wfn/relreference.h>
#include <src/integral/os/overlapbatch.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/smallints1e_london.h>
#include <src/ci/zfci/zharrison.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelReference)

using namespace std;
using namespace bagel;

shared_ptr<Reference> RelReference::project_coeff(shared_ptr<const Geometry> geomin, const bool check_geom_change) const {
  assert(check_geom_change);

  shared_ptr<Reference> out;
  const bool giao = (geomin->magnetism() || geom_->magnetism());
  if (!geomin->magnetism() && geom_->magnetism())
    throw runtime_error("Projection from GIAO to real basis set is not implemented.");

  bool moved = false;
  bool newbasis = false;
  bool newfield = (giao ? (geomin->magnetic_field() != geom_->magnetic_field()) : false);

  if (check_geom_change) {
    auto j = geomin->atoms().begin();
    for (auto& i : geom_->atoms()) {
      moved |= i->distance(*j) > 1.0e-12;
      newbasis |= i->basis() != (*j)->basis();
      ++j;
    }
  } else {
    newbasis = true;
  }

  if (moved && newbasis)
    throw runtime_error("changing geometry and basis set at the same time is not allowed");

  if (geomin->magnetism() && !geom_->magnetism())
    if (geomin->nonzero_magnetic_field() || moved || newbasis)
      throw runtime_error("The conversion from standard orbitals to GIAO requires that no simultaneous changes be made to atom positions, basis set, or magnetic field.");


  if (newbasis) {
    // 4-component wavefunction, change of basis
    // in this case we first form overlap and S^-1 matrices
    shared_ptr<ZMatrix> sinv;
    if (giao)
      sinv = make_shared<RelOverlap_London>(geomin);
    else
      sinv = make_shared<RelOverlap>(geomin);
    shared_ptr<ZMatrix> overlap = sinv->copy();
    sinv->inverse();

    const int nb = geomin->nbasis();
    const int mb = geom_->nbasis();
    ZMatrix mixed(nb*4, mb*4);

    if (!giao) {

      MixedBasis<OverlapBatch> smixed(geom_, geomin);
      MixedBasis<KineticBatch> tmixed(geom_, geomin);

      const complex<double> one(1.0);
      const complex<double> sca = one * (0.5/(c__*c__));
      mixed.copy_real_block(one,    0,    0, nb, mb, smixed);
      mixed.copy_real_block(one,   nb,   mb, nb, mb, smixed);
      mixed.copy_real_block(sca, 2*nb, 2*mb, nb, mb, tmixed);
      mixed.copy_real_block(sca, 3*nb, 3*mb, nb, mb, tmixed);
    } else {

      shared_ptr<const Geometry> relgeomin = geomin->relativistic(false);
      MixedBasis<ComplexOverlapBatch, ZMatrix> smixed(geom_, relgeomin);
      MixedBasisArray<SmallInts1e_London<ComplexOverlapBatch>, ZMatrix> smallovl(geom_, relgeomin);

      const complex<double> r2 (0.25 / (c__*c__));
      const complex<double> i2 (0.0, r2.real());
      mixed.copy_block(0,    0, nb, mb, smixed);
      mixed.copy_block(nb,  mb, nb, mb, smixed);
      mixed.add_block( r2, 2*nb, 2*mb, nb, mb, *smallovl.data(0));
      mixed.add_block( r2, 3*nb, 3*mb, nb, mb, *smallovl.data(0));
      mixed.add_block( i2, 2*nb, 2*mb, nb, mb, *smallovl.data(1));
      mixed.add_block(-i2, 3*nb, 3*mb, nb, mb, *smallovl.data(1));
      mixed.add_block( i2, 2*nb, 3*mb, nb, mb, *smallovl.data(2));
      mixed.add_block( i2, 3*nb, 2*mb, nb, mb, *smallovl.data(2));
      mixed.add_block( r2, 2*nb, 3*mb, nb, mb, *smallovl.data(3));
      mixed.add_block(-r2, 3*nb, 2*mb, nb, mb, *smallovl.data(3));
    }
    auto c = make_shared<ZMatrix>(*sinv * mixed * *relcoeff_);

    // make coefficient orthogonal
    ZMatrix unit = *c % *overlap * *c;
    unit.inverse_half();
    *c *= unit;

    auto c2 = make_shared<ZCoeff_Striped>(*c, relcoeff_->nclosed(), relcoeff_->nact(), relcoeff_->nvirt_nr(), relcoeff_->nneg());
    out = make_shared<RelReference>(geomin, c2, energy_, nneg(), nclosed(), nact(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, kramers_);

  } else if (!moved && !newfield) {
    // Special case - do nothing when converting from standard basis to GIAO without adding field
    out = make_shared<RelReference>(geomin, relcoeff_, energy_, nneg(), nclosed(), nact(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, kramers_);

  } else {

    // 4-component wavefunction, change of atom positions (or change in magnetic field)
    shared_ptr<ZMatrix> snew, sold;
    if (!giao) {
      snew = make_shared<RelOverlap>(geomin);
      sold = make_shared<RelOverlap>(geom_);
    } else {
      snew = make_shared<RelOverlap_London>(geomin);
      sold = make_shared<RelOverlap_London>(geom_);
    }
    snew->inverse_half();
    sold->sqrt();
    auto c = make_shared<ZMatrix>(*snew * *sold * *relcoeff_);

    auto c2 = make_shared<ZCoeff_Striped>(*c, relcoeff_->nclosed(), relcoeff_->nact(), relcoeff_->nvirt_nr(), relcoeff_->nneg());
    out = make_shared<RelReference>(geomin, c2, energy_, nneg(), nclosed(), nact(), nvirt()+2*(geomin->nbasis()-geom_->nbasis()), gaunt_, breit_, kramers_);
  }
  return out;
}


shared_ptr<const Kramers<2,ZRDM<1>>> RelReference::rdm1(const int ist, const int jst) const {
  assert(ciwfn_);
  ZFCI_bare fci(ciwfn_);
  return fci.rdm1(ist, jst);
}


shared_ptr<const Kramers<4,ZRDM<2>>> RelReference::rdm2(const int ist, const int jst) const {
  assert(ciwfn_);
  ZFCI_bare fci(ciwfn_);
  return fci.rdm2(ist, jst);
}


shared_ptr<const Kramers<6,ZRDM<3>>> RelReference::rdm3(const int ist, const int jst) const {
  assert(ciwfn_);
  ZFCI_bare fci(ciwfn_);
  return fci.rdm3(ist, jst);
}


shared_ptr<const Kramers<8,ZRDM<4>>> RelReference::rdm4(const int ist, const int jst) const {
  assert(ciwfn_);
  ZFCI_bare fci(ciwfn_);
  return fci.rdm4(ist, jst);
}


shared_ptr<Reference> RelReference::extract_state(const vector<int> input, const bool update_rdms) const {
  ZFCI_bare fci(ciwfn_);
  using PairType = pair<shared_ptr<const RelSpace>,shared_ptr<const RelSpace>>;
  cout << " * Extracting CI coefficients from RelReference object for the following states: ";
  for (int i = 0; i != input.size(); ++i)
    cout << input[i] << " ";
  cout << endl;

  vector<double> newenergies(input.size());
  for (int i = 0; i != input.size(); ++i)
    newenergies[i] = energy_[input[i]];

  // Construct a RelCIWfn with only CI coefficients for the desired state
  auto newciwfn = make_shared<RelCIWfn>(geom_, fci.ncore(), fci.norb(), input.size(), newenergies,
                                        ciwfn_->civectors()->extract_state(input),
                                        make_shared<PairType>(make_pair(ciwfn_->det()->first, ciwfn_->det()->second)));

  // Use extract_average_rdm(...) to get desired RDMs and prepare output
  shared_ptr<RelReference> out;
  if (update_rdms) {
    shared_ptr<RelReference> rdmref = dynamic_pointer_cast<RelReference>(extract_average_rdm(input));
    out = make_shared<RelReference>(geom_, relcoeff_, newenergies, nneg_, nclosed_, nact_, nvirt_, gaunt_, breit_,
                                    kramers_, rdmref->rdm1_av(), rdmref->rdm2_av(), newciwfn);
  } else {
    out = make_shared<RelReference>(geom_, relcoeff_, newenergies, nneg_, nclosed_, nact_, nvirt_, gaunt_, breit_,
                                    kramers_, rdm1_av(), rdm2_av(), newciwfn);
  }
  return out;
}


// TODO Cleanup or remove?  Body is mostly the same as ZHarrison::rdm12()
shared_ptr<Reference> RelReference::extract_average_rdm(const vector<int> rdm_state) const {
  ZFCI_bare fci(ciwfn_);
  if (rdm_state.size() == 0 || rdm_state.size() > nstate())
    throw runtime_error("Trying to obtain a state-averaged RDM over some invalid number of states.");

  cout << " * Extracting RDMs for ";
  cout << (rdm_state.size() > 1 ? "the average of the following states: " : "the following state: ");
  for (int i = 0; i != rdm_state.size(); ++i)
    cout << rdm_state[i] << " ";
  cout << endl;

  // for one-body RDM
  vector<shared_ptr<Kramers<2,ZRDM<1>>>> rdm1;
  vector<shared_ptr<Kramers<4,ZRDM<2>>>> rdm2;
  rdm1.resize(rdm_state.size());
  rdm2.resize(rdm_state.size());

  shared_ptr<Kramers<2,ZRDM<1>>> rdm1_av;
  shared_ptr<Kramers<4,ZRDM<2>>> rdm2_av;

  for (int index = 0; index != rdm_state.size(); ++index) {
    const int istate = rdm_state[index];

    // one body RDM
    rdm1[index] = fci.rdm1(istate, istate);

    // if off-diagonals are zero, generate a blank RDM for completeness
    if (!rdm1[index]->exist({1,0}))
      rdm1[index]->add({1,0}, rdm1[index]->at({0,0})->clone());

    // two body RDM
    rdm2[index] = fci.rdm2(istate, istate);

    // append permutation information
    rdm2[index]->emplace_perm({{0,3,2,1}},-1);
    rdm2[index]->emplace_perm({{2,3,0,1}}, 1);
    rdm2[index]->emplace_perm({{2,1,0,3}},-1);
  }

  if (rdm_state.size() > 1) {
    rdm1_av = make_shared<Kramers<2,ZRDM<1>>>();
    rdm2_av = make_shared<Kramers<4,ZRDM<2>>>();
    for (int index = 0; index != rdm_state.size(); ++index) {
      for (auto& i : *rdm1[index])
        rdm1_av->add(i.first, i.second);
      for (auto& i : *rdm2[index])
        rdm2_av->add(i.first, i.second);
    }
    for (auto& i : *rdm1_av) i.second->scale(1.0/rdm_state.size());
    for (auto& i : *rdm2_av) i.second->scale(1.0/rdm_state.size());
    rdm2_av->emplace_perm({{0,3,2,1}},-1);
    rdm2_av->emplace_perm({{2,3,0,1}}, 1);
    rdm2_av->emplace_perm({{2,1,0,3}},-1);
  } else {
    rdm1_av = rdm1.front();
    rdm2_av = rdm2.front();
  }


  auto rdm1_out = make_shared<ZMatrix>(2*fci.norb(), 2*fci.norb());
  auto rdm2_out = make_shared<ZMatrix>(4*fci.norb()*fci.norb(), 4*fci.norb()*fci.norb());
  shared_ptr<const ZRDM<1>> tmp1 = expand_kramers<1,complex<double>>(rdm1_av, fci.norb());
  shared_ptr<const ZRDM<2>> tmp2 = expand_kramers<2,complex<double>>(rdm2_av, fci.norb());
  copy_n(tmp1->data(), tmp1->size(), rdm1_out->data());
  copy_n(tmp2->data(), tmp2->size(), rdm2_out->data());

  return make_shared<RelReference>(geom_, relcoeff_, energy_, nneg_, nclosed_, nact_, nvirt_, gaunt_, breit_,
                                   kramers_, rdm1_out, rdm2_out, ciwfn_);
}
