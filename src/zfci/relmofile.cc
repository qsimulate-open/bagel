//
// BAGEL - Parallel electron correlation program.
// Filename: relmofile.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/zfci/relmofile.h>
#include <src/rel/reloverlap.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

RelMOFile::RelMOFile(const shared_ptr<const Reference> ref, const shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> co)
 : geom_(geom), ref_(dynamic_pointer_cast<const RelReference>(ref)), coeff_(co) {
  // input should be RelReference
  assert(dynamic_pointer_cast<const RelReference>(ref));

  // density fitting is assumed
  assert(geom_->df());
}


// nstart and nfence are based on the convention in Dirac calculations
void RelMOFile::init(const int nstart, const int nfence) {
  // first compute all the AO integrals in core
  nbasis_ = geom_->nbasis();
  nocc_ = (nfence - nstart)/2;
  assert((nfence - nstart) % 2 == 0);
  if (!geom_->dfs())
    geom_ = geom_->relativistic(ref_->gaunt());

  // calculates the core fock matrix
  shared_ptr<const ZMatrix> hcore = make_shared<RelHcore>(geom_);
  if (nstart != 0) {
    shared_ptr<const ZMatrix> den = coeff_->distmatrix()->form_density_rhf(nstart)->matrix();
    core_fock_ = make_shared<DFock>(geom_, hcore, coeff_->slice(0, nstart), ref_->gaunt(), ref_->breit(), /*do_grad = */false, /*robust*/ref_->breit());
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

  // then compute Kramers adapated coefficient matrices
  kramers_coeff_ = kramers(nstart, nfence);

  // calculate 1-e MO integrals
  unordered_map<bitset<2>, shared_ptr<const ZMatrix>> buf1e = compute_mo1e(kramers_coeff_);

  // calculate 2-e MO integrals
  unordered_map<bitset<4>, shared_ptr<const ZMatrix>> buf2e = compute_mo2e(kramers_coeff_);

  // compress and set mo1e_ and mo2e_
  compress_and_set(buf1e, buf2e);
}


array<shared_ptr<const ZMatrix>,2> RelMOFile::kramers(const int nstart, const int nfence) const {
  shared_ptr<const ZMatrix> coeff = coeff_->slice(nstart, nfence);
  shared_ptr<ZMatrix> reordered = coeff->clone();
  auto overlap = make_shared<RelOverlap>(geom_);

  const int noff = reordered->mdim()/2;
  const int ndim = reordered->ndim();
  const int mdim = reordered->mdim();
  const int nb = ndim / 4;
  assert(nb == nbasis_);

  if (nfence-nstart <= 0 || (nfence-nstart)%2 != 0 || ndim%4 != 0)
    throw logic_error("illegal call of RelMOFile::kramers");

  // overlap matrix
  auto sigmaz = overlap->copy();
  sigmaz->add_block(-2.0, nb, nb, nb, nb, sigmaz->get_submatrix(nb,nb,nb,nb));
  sigmaz->add_block(-2.0, nb*3, nb*3, nb, nb, sigmaz->get_submatrix(nb*3,nb*3,nb,nb));
  // just for convenience
  sigmaz->scale(-1.0);

  unique_ptr<double[]> tmp(new double[mdim]);
  int i;
  for (i = 0; i != mdim; ) {
    const double eig = ref_->eig()[nstart+i];
    int j = i+1;
    // pick up degenerate orbitals
    while (j != mdim && (fabs(ref_->eig()[nstart+j]-eig) < 1.0e-8))
      ++j;
    assert((j-i)%2 == 0);
    const int n = j-i;

    auto cnow = coeff->slice(i, j);
    auto corig = cnow->copy();
    auto s = make_shared<ZMatrix>(*cnow % *sigmaz * *cnow);
    s->diagonalize(tmp.get());
    *cnow *= *s;

    // fix the phase - making the largest large-component element in each colomn real
    for (int i = 0; i != n; ++i) {
      const int iblock = i/(n/2);
      complex<double> ele = *max_element(cnow->element_ptr(iblock*nb,i), cnow->element_ptr((iblock+1)*nb,i),
                                         [](complex<double> a, complex<double> b) { return norm(a)+1.0e-5 < norm(b); }); // favors the first one
      const complex<double> fac = norm(ele) / ele;
      transform(cnow->element_ptr(0,i), cnow->element_ptr(0,i+1), cnow->element_ptr(0,i), [&fac](complex<double> a) { return a*fac; });
    }

    // off diagonal
    const int m = n/2;
    cnow->add_block(-1.0, nb, 0, nb, m, cnow->get_submatrix(0, m, nb, m)->get_conjg());
    cnow->copy_block(0, m, nb, m, (*cnow->get_submatrix(nb, 0, nb, m)->get_conjg() * (-1.0)));
    cnow->add_block(-1.0, nb*3, 0, nb, m, cnow->get_submatrix(nb*2, m, nb, m)->get_conjg());
    cnow->copy_block(nb*2, m, nb, m, (*cnow->get_submatrix(nb*3, 0, nb, m)->get_conjg() * (-1.0)));

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

    assert(i%2 == 0);
    reordered->copy_block(0, i/2,      ndim, n/2, cnow->element_ptr(0, 0));
    reordered->copy_block(0, i/2+noff, ndim, n/2, cnow->element_ptr(0, n/2));

    i = j;
  }

  array<shared_ptr<ZMatrix>,2> out{{reordered->slice(0,noff), reordered->slice(noff, noff*2)}};

  return array<shared_ptr<const ZMatrix>,2>{{out[0], out[1]}};
}


void RelMOFile::compress_and_set(unordered_map<bitset<2>,shared_ptr<const ZMatrix>> buf1e, unordered_map<bitset<4>,shared_ptr<const ZMatrix>> buf2e) {
  mo1e_ = buf1e;

  // Harrison requires <ij|kl> = (ik|jl)
  for (auto& mat : buf2e) {
    shared_ptr<ZMatrix> tmp = mat.second->clone();
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(mat.second->data(), tmp->data(), nocc_, nocc_, nocc_, nocc_);
    bitset<4> s = mat.first;
    s[2] = mat.first[1];
    s[1] = mat.first[2];
    mo2e_.insert(make_pair(s, tmp));
  }
}


unordered_map<bitset<2>, shared_ptr<const ZMatrix>> RelJop::compute_mo1e(const array<shared_ptr<const ZMatrix>,2> coeff) {
  unordered_map<bitset<2>, shared_ptr<const ZMatrix>> out;

  for (size_t i = 0; i != 4; ++i)
    out[bitset<2>(i)] = make_shared<ZMatrix>(*coeff[i/2] % *core_fock_ * *coeff[i%2]);

  assert(out.size() == 4);
  // symmetry requirement
  assert((*out[bitset<2>("10")] - *out[bitset<2>("01")]->transpose_conjg()).rms() < 1.0e-8);
  // Kramers requirement
  assert((*out[bitset<2>("11")] - *out[bitset<2>("00")]->get_conjg()).rms() < 1.0e-8);

  return out;
}


unordered_map<bitset<4>, shared_ptr<const ZMatrix>> RelJop::compute_mo2e(const array<shared_ptr<const ZMatrix>,2> coeff) {

  auto compute = [&coeff, this](unordered_map<bitset<4>, shared_ptr<ZMatrix>>& out, const bool gaunt, const bool breit) {
    assert(!breit || gaunt);
    // (1) make DFDists
    vector<shared_ptr<const DFDist>> dfs;
    if (!gaunt) {
      dfs = geom_->dfs()->split_blocks();
      dfs.push_back(geom_->df());
    } else {
      dfs = geom_->dfsl()->split_blocks();
    }
    list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, gaunt);

    // Separate Coefficients into real and imaginary
    // correlated occupied orbitals
    array<array<shared_ptr<const Matrix>,4>,2> rocoeff;
    array<array<shared_ptr<const Matrix>,4>,2> iocoeff;
    for (int k = 0; k != 2; ++k) {
      for (int i = 0; i != 4; ++i) {
        shared_ptr<const ZMatrix> oc = coeff[k]->get_submatrix(i*nbasis_, 0, nbasis_, nocc_);
        assert(nocc_ == coeff[k]->mdim());
        rocoeff[k][i] = oc->get_real_part();
        iocoeff[k][i] = oc->get_imag_part();
      }
    }

    // (2) first-transform
    array<list<shared_ptr<RelDFHalf>>,2> half_complex;
    for (int k = 0; k != 2; ++k)
      half_complex[k] = DFock::make_half_complex(dfdists, rocoeff[k], iocoeff[k]);

    // (3) split and factorize
    array<list<shared_ptr<RelDFHalf>>,2> half_complex_exch, half_complex_exch2;
    for (size_t k = 0; k != 2; ++k) {
      for (auto& i : half_complex[k]) {
        list<shared_ptr<RelDFHalf>> tmp = i->split(/*docopy=*/false);
        half_complex_exch[k].insert(half_complex_exch[k].end(), tmp.begin(), tmp.end());
      }
      half_complex[k].clear();
      DFock::factorize(half_complex_exch[k]);
    }

    // ** save 1 external integrals (to be used in CASSCF) *** //
    if (!gaunt) {
      half_complex_coulomb_ = half_complex_exch;
    } else {
      half_complex_gaunt_ = half_complex_exch;
    }

    if (breit) {
      // TODO Not the best implementation -- one could avoid apply_J to half-transformed objects
      auto breitint = make_shared<BreitInt>(geom_);
      list<shared_ptr<Breit2Index>> breit_2index;
      for (int i = 0; i != breitint->Nblocks(); ++i) {
        breit_2index.push_back(make_shared<Breit2Index>(breitint->index(i), breitint->data(i), geom_->df()->data2()));
        if (breitint->not_diagonal(i))
          breit_2index.push_back(breit_2index.back()->cross());
      }
      for (size_t k = 0; k != 2; ++k) {
        for (auto& i : half_complex_exch[k])
          half_complex_exch2[k].push_back(i->apply_J());

        for (auto& i : half_complex_exch[k])
          for (auto& j : breit_2index)
            if (i->alpha_matches(j)) {
              half_complex_exch2[k].push_back(i->apply_J()->multiply_breit2index(j));
              DFock::factorize(half_complex_exch2[k]);
            }
      }
    }

    // (4) compute (gamma|ii)
    unordered_map<bitset<2>, shared_ptr<const RelDFFull>> full = compute_full(rocoeff, iocoeff, half_complex_exch, true);
    unordered_map<bitset<2>, shared_ptr<const RelDFFull>> full2 = !breit ? full : compute_full(rocoeff, iocoeff, half_complex_exch2, false);

    // (5) compute 4-index quantities (16 of them - we are not using symmetry... and this is a very cheap step)
    const double gscale = gaunt ? (breit ? -0.5 : -1.0) : 1.0;
    for (size_t i = 0; i != 16; ++i) {
      // we do not need (1000, 0111, 1110, 0001, 1100, 0110)
      if (i == 8 || i == 7 || i == 14 || i == 1 || i == 12 || i == 6)
        continue;

      // we will construct (1111, 1010, 1101, 0100) later
      if (i == 15 || i == 10 || i == 13 || i == 4)
        continue;

      // we compute: 0000, 0010, 1001, 0101, 0011, 1011

      const bitset<2> b2a = bitset<2>(i/4);
      const bitset<2> b2b = bitset<2>(i%4);
      const bitset<4> b4 = bitset<4>(i);
      if (!breit) {
        if (out.find(b4) == out.end()) {
          out[b4] = full.at(b2a)->form_4index(full2.at(b2b), gscale);
        } else {
          *out.at(b4) += *full.at(b2a)->form_4index(full2.at(b2b), gscale);
        }
      } else {
        // in breit cases we explicitly symmetrize the Hamiltnian
        if (out.find(b4) == out.end()) {
          out[b4] = full.at(b2a)->form_4index(full2.at(b2b), gscale*0.5);
        } else {
          *out.at(b4) += *full.at(b2a)->form_4index(full2.at(b2b), gscale*0.5);
        }
        *out.at(b4) += *full2.at(b2a)->form_4index(full.at(b2b), gscale*0.5);
      }
    }
  };

  unordered_map<bitset<4>, shared_ptr<ZMatrix>> out;
  // Dirac-Coulomb term
  compute(out, false, false);

  if (ref_->gaunt())
    compute(out, true, ref_->breit());

  // Kramers and particle symmetry
  out[bitset<4>("1111")] = out.at(bitset<4>("0000"))->get_conjg();

  out[bitset<4>("1010")] = out.at(bitset<4>("0101"))->clone();
  shared_ptr<ZMatrix> m1010 = out.at(bitset<4>("0101"))->get_conjg();
  SMITH::sort_indices<1,0,3,2,0,1,1,1>(m1010->data(), out[bitset<4>("1010")]->data(), nocc_, nocc_, nocc_, nocc_);

  out[bitset<4>("1101")] = out.at(bitset<4>("1011"))->clone();
  shared_ptr<ZMatrix> m1101 = out.at(bitset<4>("1011"))->get_conjg();
  SMITH::sort_indices<3,2,1,0,0,1,1,1>(m1101->data(), out[bitset<4>("1101")]->data(), nocc_, nocc_, nocc_, nocc_);

  out[bitset<4>("0100")] = out.at(bitset<4>("0010"))->clone();
  shared_ptr<ZMatrix> m0100 = out.at(bitset<4>("0010"))->get_conjg();
  SMITH::sort_indices<3,2,1,0,0,1,1,1>(m0100->data(), out[bitset<4>("0100")]->data(), nocc_, nocc_, nocc_, nocc_);

#if 0
  // for completeness we can compute the others too
  vector<int> target{8, 7, 14, 1, 6};
  for (auto& t : target) {
    bitset<4> tb(t);
    bitset<4> sb; sb[0] = tb[1]; sb[1] = tb[0]; sb[2] = tb[3]; sb[3] = tb[2];
    assert(out.find(tb) == out.end());
    out[tb] = out.at(sb)->clone();
    SMITH::sort_indices<1,0,3,2,0,1,1,1>(out.at(sb)->data(), out.at(tb)->data(), nocc_, nocc_, nocc_, nocc_);
    transform(out.at(tb)->data(), out.at(tb)->data()+nocc_*nocc_*nocc_*nocc_, out.at(tb)->data(), [](complex<double> a) { return conj(a); });
  }
  out[bitset<4>("1100")] = out.at(bitset<4>("0011"))->transpose();
#endif

  return unordered_map<bitset<4>, shared_ptr<const ZMatrix>>(out.begin(), out.end());
}
