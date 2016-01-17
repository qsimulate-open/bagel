//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relmofile.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
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

#include <iostream>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/util/prim_op.h>
#include <src/ci/zfci/relmofile.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/scf/dhf/dfock.h>

using namespace std;
using namespace bagel;

RelMOFile::RelMOFile(const shared_ptr<const Geometry> geom, shared_ptr<const RelCoeff_Block> co, const int charge, const bool gaunt, const bool breit, const bool tsymm)
 : charge_(charge), geom_(geom), coeff_(co), gaunt_(gaunt), breit_(breit), tsymm_(tsymm) {
  // density fitting is assumed
  assert(geom_->df());
}


// nstart and nfence are based on the convention in Dirac calculations
void RelMOFile::init(const int nstart, const int nfence, const bool restricted) {
  // first compute all the AO integrals in core
  nbasis_ = geom_->nbasis();
  nocc_ = (nfence - nstart)/2;
  assert((nfence - nstart) % 2 == 0);
  assert(geom_->dfs());

  // calculates the core fock matrix
  shared_ptr<const ZMatrix> hcore;
  if (!geom_->magnetism())
    hcore = make_shared<RelHcore>(geom_);
  else
    hcore = make_shared<RelHcore_London>(geom_);

  if (nstart != 0) {
    shared_ptr<const ZMatrix> den = coeff_->distmatrix()->form_density_rhf(nstart)->matrix();
    core_fock_ = make_shared<DFock>(geom_, hcore, coeff_->slice_copy(0, nstart), gaunt_, breit_, /*do_grad = */false, /*robust*/breit_);
    const complex<double> prod = (*den * (*hcore+*core_fock_)).trace();
    if (fabs(prod.imag()) > 1.0e-12) {
      stringstream ss; ss << "imaginary part of energy is nonzero!! Perhaps Fock is not Hermite for some reasons " << setprecision(10) << prod.imag();
      cout << ss.str() << endl;
    }
    core_energy_ = 0.5*prod.real();
  } else {
    core_fock_ = hcore;
    core_energy_ = 0.0;
  }

  kramers_coeff_ = coeff_->kramers_active();

  // calculate 1-e MO integrals
  shared_ptr<Kramers<2,ZMatrix>> buf1e = compute_mo1e(kramers_coeff_);

  // calculate 2-e MO integrals
  shared_ptr<Kramers<4,ZMatrix>> buf2e = compute_mo2e(kramers_coeff_);

  // compress and set mo1e_ and mo2e_
  compress_and_set(buf1e, buf2e);
}



// this is a static function!
shared_ptr<Kramers<1,ZMatrix>> RelMOFile::kramers(shared_ptr<const ZMatrix> coeff, shared_ptr<const ZMatrix> overlap, shared_ptr<const ZMatrix> hcore) {
  const int noff = coeff->mdim()/2;
  const int ndim = coeff->ndim();
  const int mdim = coeff->mdim();
  const int nb = ndim / 4;
  unique_ptr<complex<double>[]> eig = (*coeff % *hcore * *coeff).diag();

  auto out = make_shared<Kramers<1,ZMatrix>>();
  out->emplace(0, make_shared<ZMatrix>(ndim, noff));
  out->emplace(1, make_shared<ZMatrix>(ndim, noff));

  if (ndim%2 != 0 || ndim%4 != 0)
    throw logic_error("illegal call of RelMOFile::kramers");

  // overlap matrix
  auto sigmaz = overlap->copy();
  sigmaz->add_block(-2.0, nb, nb, nb, nb, sigmaz->get_submatrix(nb,nb,nb,nb));
  sigmaz->add_block(-2.0, nb*3, nb*3, nb, nb, sigmaz->get_submatrix(nb*3,nb*3,nb,nb));
  // just for convenience
  sigmaz->scale(-1.0);

  VectorB tmp(mdim);

  list<int> done;
  for (int i = 0; i != mdim; ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    list<int> current{i};
    const double e = eig[i].real();

    for (int j = i+1; j < mdim; ++j) {
      if (fabs(eig[j].real()-e)/fabs(e) < 1.0e-8)
        current.push_back(j);
    }
    const int n = current.size();
    if (n%2 != 0) throw runtime_error("orbitals are not kramers paired");

    auto cnow = make_shared<ZMatrix>(ndim, n);
    int j = 0;
    for (auto& i : current)
      cnow->copy_block(0, j++, ndim, 1, coeff->element_ptr(0,i));

    auto corig = cnow->copy();
    auto s = make_shared<ZMatrix>(*cnow % *sigmaz * *cnow);
    s->diagonalize(tmp);
    *cnow *= *s;

    // fix the phase - making the largest large-component element in each column real
#if 1
    for (int i = 0; i != n; ++i) {
      const int iblock = i/(n/2);
      complex<double> ele = *max_element(cnow->element_ptr(iblock*nb,i), cnow->element_ptr((iblock+1)*nb,i),
                                         [](complex<double> a, complex<double> b) { return norm(a)+1.0e-5 < norm(b); }); // favors the first one
      const complex<double> fac = norm(ele) / ele*complex<double>(1.0,1.0);
      for_each(cnow->element_ptr(0,i), cnow->element_ptr(0,i+1), [&fac](complex<double>& a) { a *= fac; });
    }
#endif

    // off diagonal
    const int m = n/2;
    cnow->add_block(-1.0, nb, 0, nb, m, cnow->get_submatrix(0, m, nb, m)->get_conjg());
    cnow->copy_block(0, m, nb, m, *cnow->get_submatrix(nb, 0, nb, m)->get_conjg() * (-1.0));
    cnow->add_block(-1.0, nb*3, 0, nb, m, cnow->get_submatrix(nb*2, m, nb, m)->get_conjg());
    cnow->copy_block(nb*2, m, nb, m, *cnow->get_submatrix(nb*3, 0, nb, m)->get_conjg() * (-1.0));

    // diagonal
    cnow->add_block(1.0, 0, 0, nb, m, cnow->get_submatrix(nb, m, nb, m)->get_conjg());
    cnow->copy_block(nb, m, nb, m, cnow->get_submatrix(0, 0, nb, m)->get_conjg());
    cnow->add_block(1.0, nb*2, 0, nb, m, cnow->get_submatrix(nb*3, m, nb, m)->get_conjg());
    cnow->copy_block(nb*3, m, nb, m, cnow->get_submatrix(nb*2, 0, nb, m)->get_conjg());

    auto diag = (*cnow % *overlap * *cnow).diag();
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndim; ++j)
        cnow->element(j,i) /= sqrt(diag[i].real());

    ZMatrix unit = *corig % *overlap * *cnow;
    unit.purify_unitary();
    *cnow = *corig * unit;

    const int d = done.size();
    assert(d % 2 == 0);
    out->at(0)->copy_block(0, d/2, ndim, n/2, cnow->element_ptr(0, 0));
    out->at(1)->copy_block(0, d/2, ndim, n/2, cnow->element_ptr(0, n/2));

    done.insert(done.end(), current.begin(), current.end());
  }
  return out;
}


void RelMOFile::compress_and_set(shared_ptr<Kramers<2,ZMatrix>> buf1e, shared_ptr<Kramers<4,ZMatrix>> buf2e) {
  mo1e_ = buf1e;
  mo2e_ = make_shared<Kramers<4,ZMatrix>>();

  // Harrison requires <ij|kl> = (ik|jl)
  for (auto& mat : *buf2e) {
    shared_ptr<ZMatrix> tmp = mat.second->clone();
    sort_indices<0,2,1,3,0,1,1,1>(mat.second->data(), tmp->data(), nocc_, nocc_, nocc_, nocc_);
    bitset<4> s = mat.first.tag();
    s[2] = mat.first.tag()[1];
    s[1] = mat.first.tag()[2];
    mo2e_->emplace(s, tmp);
  }
}


shared_ptr<Kramers<2,ZMatrix>> RelJop::compute_mo1e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {
  auto out = make_shared<Kramers<2,ZMatrix>>();
  const int n = tsymm_ ? 3 : 4;
  for (size_t i = 0; i != n; ++i)
    out->emplace(i, make_shared<ZMatrix>(*coeff->at(i/2) % *core_fock_ * *coeff->at(i%2)));
  if (tsymm_)
    out->emplace({1,1}, out->at({0,0})->get_conjg());

  assert(out->size() == 4);
  // symmetry requirement
  assert((*out->at({1,0}) - *out->at({0,1})->transpose_conjg()).rms() < 1.0e-8);
  // Kramers requirement
  assert((*make_shared<ZMatrix>(*coeff->at(1) % *core_fock_ * *coeff->at(1)) - *out->at({0,0})->get_conjg()).rms() < 1.0e-8 || !tsymm_);

  return out;
}


tuple<list<shared_ptr<RelDFHalf>>,list<shared_ptr<RelDFHalf>>>
  RelMOFile::compute_half(shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> coeff, const bool gaunt, const bool breit) {

  assert(!breit || gaunt);
  // (1) make DFDists
  vector<shared_ptr<const DFDist>> dfs;
  if (!gaunt) {
    dfs = geom->dfs()->split_blocks();
    dfs.push_back(geom->df());
  } else {
    dfs = geom->dfsl()->split_blocks();
  }
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, gaunt);

  // (2) first-transform
  list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, coeff);

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_exch, half_complex_exch2;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(/*docopy=*/false);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }
  half_complex.clear();
  DFock::factorize(half_complex_exch);

  if (breit) {
    // TODO Not the best implementation -- one could avoid apply_J to half-transformed objects
    auto breitint = make_shared<BreitInt>(geom);
    list<shared_ptr<Breit2Index>> breit_2index;
    for (int i = 0; i != breitint->Nblocks(); ++i) {
      breit_2index.push_back(make_shared<Breit2Index>(breitint->index(i), breitint->data(i), geom->df()->data2()));
      if (breitint->not_diagonal(i))
        breit_2index.push_back(breit_2index.back()->cross());
    }
    for (auto& i : half_complex_exch)
      half_complex_exch2.push_back(i->apply_J());

    for (auto& i : half_complex_exch)
      for (auto& j : breit_2index)
        if (i->alpha_matches(j)) {
          half_complex_exch2.push_back(i->apply_J()->multiply_breit2index(j));
          DFock::factorize(half_complex_exch2);
        }
  }
  return make_tuple(half_complex_exch, half_complex_exch2);
}


shared_ptr<Kramers<4,ZMatrix>> RelJop::compute_mo2e(shared_ptr<const Kramers<1,ZMatrix>> coeff) {

  auto compute = [&coeff, this](shared_ptr<Kramers<4,ZMatrix>> out, const bool gaunt, const bool breit) {
    array<list<shared_ptr<RelDFHalf>>,2> half_complex_exch, half_complex_exch2;
    for (int k = 0; k != 2; ++k)
      tie(half_complex_exch[k], half_complex_exch2[k]) = compute_half(geom_, coeff->at(k), gaunt, breit);

    if (!gaunt)
      half_complex_coulomb_ = half_complex_exch;
    else
      half_complex_gaunt_ = half_complex_exch;

    // (4) compute (gamma|ii)
    auto full = make_shared<Kramers<2,ListRelDFFull>>();
    full->emplace({0,0}, compute_full(coeff->at(0), half_complex_exch[0], true));
    full->emplace({1,0}, compute_full(coeff->at(0), half_complex_exch[1], true));
    full->emplace({0,1}, compute_full(coeff->at(1), half_complex_exch[0], true));
    full->emplace({1,1}, compute_full(coeff->at(1), half_complex_exch[1], true));

    auto full2 = !breit ? full : make_shared<Kramers<2,ListRelDFFull>>();
    if (breit) {
      full2->emplace({0,0}, compute_full(coeff->at(0), half_complex_exch2[0], false));
      full2->emplace({1,0}, compute_full(coeff->at(0), half_complex_exch2[1], false));
      full2->emplace({0,1}, compute_full(coeff->at(1), half_complex_exch2[0], false));
      full2->emplace({1,1}, compute_full(coeff->at(1), half_complex_exch2[1], false));
    }

    // (5) compute 4-index quantities (16 of them - we are not using symmetry... and this is a very cheap step)
    const double gscale = gaunt ? (breit ? -0.25 : -1.0) : 1.0;
    for (size_t i = 0; i != 16; ++i) {
      // we do not need (1000, 0111, 1110, 0001, 1100, 0110)
      if (i == 8 || i == 7 || i == 14 || i == 1 || i == 12 || i == 6)
        continue;

      // we will construct (1111, 1010, 1101, 0100) later
      if ((tsymm_ && i == 15) || i == 10 || i == 13 || i == 4)
        continue;

      // we compute: 0000, 0010, 1001, 0101, 0011, 1011

      // TODO : put in if statement for apply_J if nact*nact is much smaller than number of MPI processes
      const int b2a = i/4;
      const int b2b = i%4;
      out->add(i, full->at(b2a)->form_4index(full2->at(b2b), gscale));

      // in breit cases we explicitly symmetrize the Hamiltnian (hence the prefactor 0.5 above)
      if (breit)
        *out->at(i) += *full2->at(b2a)->form_4index(full->at(b2b), gscale);
    }
  };

  auto out = make_shared<Kramers<4,ZMatrix>>();
  // Dirac-Coulomb term
  compute(out, false, false);

  if (gaunt_)
    compute(out, true, breit_);

  // Kramers and particle symmetry
  if (tsymm_)
    (*out)[{1,1,1,1}] = out->at({0,0,0,0})->get_conjg();

  (*out)[{1,0,1,0}] = out->at({0,1,0,1})->clone();
  shared_ptr<ZMatrix> m1010 = out->at({0,1,0,1})->get_conjg();
  sort_indices<1,0,3,2,0,1,1,1>(m1010->data(), out->at({1,0,1,0})->data(), nocc_, nocc_, nocc_, nocc_);

  (*out)[{1,1,0,1}] = out->at({1,0,1,1})->clone();
  shared_ptr<ZMatrix> m1101 = out->at({1,0,1,1})->get_conjg();
  sort_indices<3,2,1,0,0,1,1,1>(m1101->data(), out->at({1,1,0,1})->data(), nocc_, nocc_, nocc_, nocc_);

  (*out)[{0,1,0,0}] = out->at({0,0,1,0})->clone();
  shared_ptr<ZMatrix> m0100 = out->at({0,0,1,0})->get_conjg();
  sort_indices<3,2,1,0,0,1,1,1>(m0100->data(), out->at({0,1,0,0})->data(), nocc_, nocc_, nocc_, nocc_);

#if 0
  // for completeness we can compute the others too
  vector<int> target{8, 7, 14, 1, 6};
  for (auto& t : target) {
    bitset<4> tb(t);
    bitset<4> sb; sb[0] = tb[1]; sb[1] = tb[0]; sb[2] = tb[3]; sb[3] = tb[2];
    assert(out.find(tb) == out.end());
    out[tb] = out.at(sb)->clone();
    sort_indices<1,0,3,2,0,1,1,1>(out.at(sb)->data(), out.at(tb)->data(), nocc_, nocc_, nocc_, nocc_);
    transform(out.at(tb)->data(), out.at(tb)->data()+nocc_*nocc_*nocc_*nocc_, out.at(tb)->data(), [](complex<double> a) { return conj(a); });
  }
  out[bitset<4>("1100")] = out.at(bitset<4>("0011"))->transpose();
#endif

  return out;
}


shared_ptr<ListRelDFFull> RelMOFile::compute_full(shared_ptr<const ZMatrix> coeff, list<shared_ptr<RelDFHalf>> half, const bool appj) {
  // TODO remove once DFDistT class is fixed
  const bool transform_with_full = !(half.front()->nocc()*coeff->mdim() <= mpi__->size());
  if (!transform_with_full && appj) {
    for (auto& i : half)
      i = i->apply_J();
  }

  list<shared_ptr<RelDFFull>> dffull;
  for (auto& i : half)
    dffull.push_back(make_shared<RelDFFull>(i, coeff));
  DFock::factorize(dffull);

  assert(dffull.size() == 1 || dffull.size() == 3);

  for (auto& i : dffull)
    i->scale(i->fac());

  if (transform_with_full && appj)
    for (auto& i : dffull)
      i = i->apply_J();
  return make_shared<ListRelDFFull>(dffull);
}


