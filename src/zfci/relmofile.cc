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

RelMOFile::RelMOFile(const shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> co, const bool gaunt, const bool breit)
 : geom_(geom), coeff_(co), gaunt_(gaunt), breit_(breit) {
  // density fitting is assumed
  assert(geom_->df());
}


// nstart and nfence are based on the convention in Dirac calculations
void RelMOFile::init(const int nstart, const int nfence, const bool restricted) {
  // first compute all the AO integrals in core
  nbasis_ = geom_->nbasis();
  nocc_ = (nfence - nstart)/2;
  assert((nfence - nstart) % 2 == 0);
  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  // calculates the core fock matrix
  shared_ptr<const ZMatrix> hcore = make_shared<RelHcore>(geom_);
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

  // then compute Kramers adapated coefficient matrices
  auto overlap = make_shared<RelOverlap>(geom_);
  if (!restricted) {
    kramers_coeff_ = kramers_zquat(nstart, nfence, coeff_->slice_copy(nstart, nfence), overlap, hcore);
  } else {
#if 1
    array<shared_ptr<const ZMatrix>,2> tmp;
    tmp[0] = make_shared<const ZMatrix>(*coeff_->slice_copy(nstart, nstart+(nfence-nstart)/2));
    tmp[1] = make_shared<const ZMatrix>(*coeff_->slice_copy(nstart+(nfence-nstart)/2, nfence));
    kramers_coeff_ = tmp;
#else
    kramers_coeff_ = kramers(coeff_->slice_copy(nstart, nfence), overlap, core_fock_);
#endif
  }

  // calculate 1-e MO integrals
  unordered_map<bitset<2>, shared_ptr<const ZMatrix>> buf1e = compute_mo1e(kramers_coeff_);

  // calculate 2-e MO integrals
  unordered_map<bitset<4>, shared_ptr<const ZMatrix>> buf2e = compute_mo2e(kramers_coeff_);

  // compress and set mo1e_ and mo2e_
  compress_and_set(buf1e, buf2e);
}



// this is a static function!
array<shared_ptr<const ZMatrix>,2> RelMOFile::kramers(shared_ptr<const ZMatrix> coeff, shared_ptr<const ZMatrix> overlap, shared_ptr<const ZMatrix> hcore) {
  const int noff = coeff->mdim()/2;
  const int ndim = coeff->ndim();
  const int mdim = coeff->mdim();
  const int nb = ndim / 4;
  unique_ptr<complex<double>[]> eig = (*coeff % *hcore * *coeff).diag();

  array<shared_ptr<ZMatrix>,2> out{{make_shared<ZMatrix>(ndim, noff), make_shared<ZMatrix>(ndim, noff)}};

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
    out[0]->copy_block(0, d/2, ndim, n/2, cnow->element_ptr(0, 0));
    out[1]->copy_block(0, d/2, ndim, n/2, cnow->element_ptr(0, n/2));

    done.insert(done.end(), current.begin(), current.end());
  }

  return array<shared_ptr<const ZMatrix>,2>{{out[0], out[1]}};
}


// modified from init_kramers_coeff in zcasscf_coeff.cc
array<shared_ptr<const ZMatrix>,2> RelMOFile::kramers_zquat(const int nstart, const int nfence, shared_ptr<const ZMatrix> coeff, shared_ptr<const ZMatrix> overlap, shared_ptr<const ZMatrix> hcore) {
  assert(coeff->mdim() > 2 && coeff->mdim()%2 == 0); // zquatev has a bug for 2x2 case since there are no super-offdiagonals in a 2x2 and tridiagonalization is probably not possible
  assert(nstart < nfence);
  const int ndim    = coeff->ndim();
  const int nb      = ndim/4;
  const int nclosed = nstart/2;
  const int nact    = (nfence - nstart)/2;
  const int nvnr    = nb - nact - nclosed;
  const int nvirt   = nvnr + nb;
  if(ndim%2 != 0 || ndim % 4 != 0)
    throw logic_error("illegal call of RelMOFile::kramers_zquat");
  // local function to transform from kramers to quaternion format
  auto quaternion = [](shared_ptr<ZMatrix> o) {
    shared_ptr<ZMatrix> scratch = o->clone();
    const int n = o->ndim()/4;
    const int m = o->mdim()/4;
    map<int, int> trans {{0,0}, {1,2}, {2,1}, {3,3}};
    for (auto& i : trans)
      for (auto& j : trans)
        scratch->copy_block(i.first*n, j.first*m, n, m, o->get_submatrix(i.second*n, j.second*m, n, m));
    *o = *scratch;
  };

  const int nocc    = nact + nclosed;

  shared_ptr<ZMatrix> focktmp;
  shared_ptr<ZMatrix> ctmp;
  {
    const int norb = geom_->nele();
    assert(norb <= coeff_->mdim());
    focktmp = make_shared<DFock>(geom_, hcore, coeff_->slice_copy(0, norb), gaunt_, breit_, /*store_half*/false, /*robust*/false);
  }
  quaternion(focktmp);

  shared_ptr<ZMatrix> s12 = overlap->tildex(1.0e-10);
  quaternion(s12);

  auto fock_tilde = make_shared<ZMatrix>(*s12 % (*focktmp) * *s12);

  // quaternion diagonalization
  {
    unique_ptr<double[]> eig(new double[fock_tilde->ndim()]);
    zquatev_(fock_tilde->ndim(), fock_tilde->data(), eig.get());
  }
  // re-order to kramers format and move negative energy states to virtual space
  ctmp = make_shared<ZMatrix>(*s12 * *fock_tilde);
  { // rows: {L+, S+, L-, S-} -> {L+, L-, S+, S-}
    assert(ctmp->ndim() % 4 == 0);
    const int n = ctmp->ndim()/4;
    const int m = ctmp->mdim();
    shared_ptr<ZMatrix> scratch = ctmp->get_submatrix(n, 0, n*2, m);
    ctmp->copy_block(n,   0, n, m, scratch->get_submatrix(n, 0, n, m));
    ctmp->copy_block(n*2, 0, n, m, scratch->get_submatrix(0, 0, n, m));
  }
  // move_positronic_orbitals
  {
    auto move_one = [this, &ctmp](const int offset, const int block1, const int block2) {
      shared_ptr<ZMatrix> scratch = make_shared<ZMatrix>(ctmp->ndim(), block1+block2);
      scratch->copy_block(0,      0, ctmp->ndim(), block2, ctmp->slice(offset+block1, offset+block1+block2));
      scratch->copy_block(0, block2, ctmp->ndim(), block1, ctmp->slice(offset,        offset+block1));
      ctmp->copy_block(0, offset, ctmp->ndim(), block1+block2, scratch);
    };
    const int nneg2 = ctmp->mdim()/4;
    move_one(           0, nneg2, nocc+nvnr);
    move_one(nocc + nvirt, nneg2, nocc+nvnr);
    { // DEBUG : TR Symmetric coeff ?
      shared_ptr<ZMatrix> ctmp2 = ctmp->copy();
      { // rows: {L+, L-, S+, S-} -> {L+, S+, L-, S-}
        assert(ctmp2->ndim() % 4 == 0);
        const int n = ctmp2->ndim()/4;
        const int m = ctmp2->mdim();
        shared_ptr<ZMatrix> scratch = ctmp2->get_submatrix(n, 0, n*2, m);
        ctmp2->copy_block(n,   0, n, m, scratch->get_submatrix(n, 0, n, m));
        ctmp2->copy_block(n*2, 0, n, m, scratch->get_submatrix(0, 0, n, m));
        cout << setprecision(8) << scientific << " coeff diag tr symm rms = " <<
          (*ctmp2->get_submatrix(0,0,n*2,m/2) - *ctmp2->get_submatrix(n*2, ctmp2->mdim()/2, n*2, m/2)->get_conjg()).rms() << endl;
        cout << setprecision(8) << scientific << " coeff off-diag tr symm rms = " <<
          (*ctmp2->get_submatrix(n*2,0,n*2,m/2) + *ctmp2->get_submatrix(0, ctmp2->mdim()/2, n*2, m/2)->get_conjg()).rms() << endl;
      }
    }

    auto tmp = coeff->clone();
    tmp->copy_block(0, 0, ctmp->ndim(), coeff->mdim()/2, ctmp->slice(nstart/2, nstart/2+coeff->mdim()/2));
    tmp->copy_block(0, coeff->mdim()/2, ctmp->ndim(), coeff->mdim()/2, ctmp->slice(ctmp->mdim()/2 + nstart/2, ctmp->mdim()/2 + nstart/2+coeff->mdim()/2));
    ctmp = tmp;
  }

  array<shared_ptr<const ZMatrix>,2> tmp;
  tmp[0] = make_shared<const ZMatrix>(*ctmp->slice_copy(0, ctmp->mdim()/2));
  tmp[1] = make_shared<const ZMatrix>(*ctmp->slice_copy(ctmp->mdim()/2, ctmp->mdim()));

  return tmp;
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
    mo2e_.emplace(s, tmp);
  }
}


unordered_map<bitset<2>, shared_ptr<const ZMatrix>> RelJop::compute_mo1e(const array<shared_ptr<const ZMatrix>,2> coeff) {
  unordered_map<bitset<2>, shared_ptr<const ZMatrix>> out;

  for (size_t i = 0; i != 4; ++i)
    out[bitset<2>(i)] = make_shared<ZMatrix>(*coeff[i/2] % *core_fock_ * *coeff[i%2]);
  out[bitset<2>("11")] = out[bitset<2>("00")]->get_conjg();

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

  if (gaunt_)
    compute(out, true, breit_);

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


void RelMOFile::update_kramers_coeff(shared_ptr<ZMatrix> coeff) {
  // transform kramers coeff to natural orbital rep
  assert(coeff->mdim() == coeff->ndim());
  assert(coeff->mdim() == nocc_*2);
  int ndim = kramers_coeff(0)->ndim();
  int n    = coeff->ndim();
  auto kcoefftmp  = make_shared<ZMatrix>(ndim, coeff->mdim());
  auto kcoefftmp2 = kcoefftmp->clone();
  kcoefftmp->copy_block(0,     0, ndim, nocc_, kramers_coeff(0)->data());
  kcoefftmp->copy_block(0, nocc_, ndim, nocc_, kramers_coeff(1)->data());
  zgemm3m_("N", "N", ndim, n, n, 1.0, kcoefftmp->data(), ndim, coeff->get_conjg()->data(), n, 0.0, kcoefftmp2->data(), ndim);
  array<shared_ptr<const ZMatrix>,2> kctmp;
  kctmp[0] = make_shared<const ZMatrix>(*kcoefftmp2->slice_copy(0, nocc_));
  kctmp[1] = make_shared<const ZMatrix>(*kcoefftmp2->slice_copy(nocc_, nocc_*2));
  kramers_coeff_ = kctmp;
}
