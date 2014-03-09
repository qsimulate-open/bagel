//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <set>
#include <src/smith/prim_op.h>
#include <src/nevpt2/nevpt2.h>
#include <src/df/dfdistt.h>
#include <src/casscf/casbfgs.h>
#include <src/casscf/qvec.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

// TODO
const static int istate = 0;

NEVPT2::NEVPT2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  casscf_ = make_shared<CASBFGS>(input, g, ref);
  casscf_->compute();
  ref_ = casscf_->conv_to_ref();

  cout << endl << "  === DF-NEVPT2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));
  norm_thresh_ = idata_->get<double>("norm_thresh", 1.0e-13);

}


void NEVPT2::compute() {

  const size_t nclosed = ref_->nclosed() - ncore_;
  const size_t nact = ref_->nact();
  const size_t nvirt = ref_->nvirt();

  // helper functions
  auto id2 = [&nact](                          const int k, const int l) { return        (       (k+nact*l)); };
  auto id3 = [&nact](             const int j, const int k, const int l) { return        (j+nact*(k+nact*l)); };
  auto id4 = [&nact](const int i, const int j, const int k, const int l) { return i+nact*(j+nact*(k+nact*l)); };

  if (nclosed+nact < 1) throw runtime_error("no correlated electrons");
  if (nvirt < 1)        throw runtime_error("no virtuals orbitals");

  // coefficients -- will be updated later
  shared_ptr<Matrix> ccoeff = ref_->coeff()->slice(ncore_, ncore_+nclosed);
  shared_ptr<Matrix> acoeff = ref_->coeff()->slice(ncore_+nclosed, ncore_+nclosed+nact);
  shared_ptr<Matrix> vcoeff = ref_->coeff()->slice(ncore_+nclosed+nact, ncore_+nclosed+nact+nvirt);
  // rdm 1
  shared_ptr<const Matrix> rdm1 = casscf_->fci()->rdm1(istate)->rdm1_mat(/*nclosed*/0);
  shared_ptr<Matrix> unit = rdm1->clone(); unit->unit();
  shared_ptr<const Matrix> hrdm1 = make_shared<Matrix>(*unit*2.0 - *rdm1);
  // rdm 2
  shared_ptr<Matrix> rdm2 = make_shared<Matrix>(nact*nact, nact*nact);
  {
    shared_ptr<const RDM<2>> r2 = ref_->rdm2(istate);
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(r2->data(), rdm2->data(), nact, nact, nact, nact);
  }
  shared_ptr<Matrix> hrdm2 = rdm2->copy();
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nact; ++j) {
      for (int k = 0; k != nact; ++k) {
        hrdm2->element(j+nact*k, i+nact*k) += 2.0 * hrdm1->element(j,i);
        hrdm2->element(k+nact*j, i+nact*k) -= hrdm1->element(j,i);
        hrdm2->element(j+nact*k, k+nact*i) +=  rdm1->element(i,j);
        hrdm2->element(k+nact*j, k+nact*i) -= 2.0 *  rdm1->element(i,j);
      }
    }
  }
  // helper matrix for higher-order hole RDMs
  shared_ptr<Matrix> srdm2 = hrdm2->clone(); // S(a,b,c,d) = <0|a+p bp cq d+q|0>
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          srdm2->element(l+nact*k,j+nact*i) = -rdm2->element(l+nact*i,k+nact*j) + (i == j ? 2.0*rdm1->element(l,k) : 0.0) - (i == k ? rdm1->element(l,j) : 0.0);

  // rdm 3 and 4
  shared_ptr<Matrix> rdm3 = make_shared<Matrix>(nact*nact*nact, nact*nact*nact);
  shared_ptr<Matrix> rdm4 = make_shared<Matrix>(nact*nact*nact*nact, nact*nact*nact*nact);
  {
    shared_ptr<const RDM<3>> r3;
    shared_ptr<const RDM<4>> r4;
    tie(r3, r4) = casscf_->fci()->compute_rdm34(istate);
    SMITH::sort_indices<0,2,4,  1,3,5,  0,1,1,1>(r3->data(), rdm3->data(), nact, nact, nact, nact, nact, nact);
    SMITH::sort_indices<0,2,4,6,1,3,5,7,0,1,1,1>(r4->data(), rdm4->data(), nact, nact, nact, nact, nact, nact, nact, nact);
  }
  shared_ptr<Matrix> hrdm3 = make_shared<Matrix>(*rdm3 * (-1.0));
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          for (int m = 0; m != nact; ++m) {
            hrdm3->element(id3(l,k,m),id3(j,i,m)) += 2.0*hrdm2->element(l+nact*k,j+nact*i);
            hrdm3->element(id3(l,m,k),id3(j,i,m)) -=     hrdm2->element(l+nact*k,j+nact*i);
            hrdm3->element(id3(m,l,k),id3(j,i,m)) -=     hrdm2->element(l+nact*k,i+nact*j);
            hrdm3->element(id3(l,k,m),id3(j,m,i)) +=     srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(l,m,k),id3(j,m,i)) -= 2.0*srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(m,l,k),id3(j,m,i)) +=     srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(l,k,m),id3(m,j,i)) -=      rdm2->element(i+nact*j,l+nact*k);
            hrdm3->element(id3(l,m,k),id3(m,j,i)) -=      rdm2->element(i+nact*j,k+nact*l);
            hrdm3->element(id3(m,l,k),id3(m,j,i)) += 2.0* rdm2->element(i+nact*j,k+nact*l);
          }
  // <a+ a b+ b> and <a+ a b+ b c+ c>
  shared_ptr<Matrix> ardm2 = rdm2->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        for (int l = 0; l != nact; ++l)
          ardm2->element(l+nact*k,j+nact*i) += rdm2->element(l+nact*j,k+nact*i);
        ardm2->element(k+nact*j,j+nact*i) += rdm1->element(k,i);
      }
  shared_ptr<Matrix> ardm3 = rdm3->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          for (int m = 0; m != nact; ++m) {
            for (int n = 0; n != nact; ++n)
              ardm3->element(id3(n,m,l),id3(k,j,i)) += rdm3->element(id3(n,l,j),id3(m,k,i));
            ardm3->element(id3(m,l,l),id3(k,j,i)) += ardm2->element(m+nact*k,j+nact*i);
            ardm3->element(id3(m,l,k),id3(j,j,i)) += rdm2->element(m+nact*k,l+nact*i);
            ardm3->element(id3(m,l,k),id3(j,l,i)) += rdm2->element(m+nact*k,i+nact*j);
          }
  shared_ptr<Matrix> ardm4 = rdm4->clone();
  for (int h = 0; h != nact; ++h)
    for (int g = 0; g != nact; ++g)
      for (int f = 0; f != nact; ++f)
        for (int e = 0; e != nact; ++e)
          for (int d = 0; d != nact; ++d)
            for (int c = 0; c != nact; ++c)
              for (int b = 0; b != nact; ++b)
                for (int a = 0; a != nact; ++a) {
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == c ? 1.0 : 0.0) * ardm3->element(id3(a,d,e),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (d == e && b == c ? 1.0 : 0.0) * ardm2->element(id2(a,f),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == e ? 1.0 : 0.0) * ardm3->element(id3(a,b,c),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (b == e && c == f ? 1.0 : 0.0) * ardm2->element(id2(a,d),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == e ? 1.0 : 0.0) * ardm3->element(id3(a,f,c),id3(d,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (f == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(b,d,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(b,h,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(h,d,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += rdm4->element(id4(a,c,e,g),id4(b,d,f,h));
                }

  // Hcore
  shared_ptr<const Matrix> hcore = make_shared<Hcore>(geom_);

  // make canonical orbitals in closed and virtual subspaces
  vector<double> veig(nvirt);
  vector<double> oeig(nclosed);
  shared_ptr<Matrix> coeffall;
  // fock
  shared_ptr<const Matrix> fockact,   fock;
  // core fock
  shared_ptr<const Matrix> fockact_c, fock_c;
  // heff_p in active
  shared_ptr<const Matrix> fockact_p, fock_p;
  // heff_h in active (active treated as closed)
  shared_ptr<const Matrix> fockact_h, fock_h;
  {
    // * core Fock operator
    shared_ptr<const Matrix> ofockao = nclosed+ncore_ ? make_shared<const Fock<1>>(geom_, hcore, nullptr, ref_->coeff()->slice(0, ncore_+nclosed), /*store*/false, /*rhf*/true) : hcore;
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> rdm1mat = rdm1->copy();
    rdm1mat->sqrt();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1mat);
    auto fockao = make_shared<Fock<1>>(geom_, ofockao, nullptr, acoeffw, /*store*/false, /*rhf*/true);
    // MO Fock
    {
      Matrix omofock(*ccoeff % *fockao * *ccoeff);
      omofock.diagonalize(oeig.data());
      *ccoeff *= omofock;
    } {
      Matrix vmofock(*vcoeff % *fockao * *vcoeff);
      vmofock.diagonalize(veig.data());
      *vcoeff *= vmofock;
    }
    coeffall = make_shared<Matrix>(ccoeff->ndim(), nclosed+nact+nvirt);
    coeffall->copy_block(0, 0           , ccoeff->ndim(), nclosed, ccoeff);
    coeffall->copy_block(0, nclosed     , acoeff->ndim(), nact   , acoeff);
    coeffall->copy_block(0, nclosed+nact, vcoeff->ndim(), nvirt  , vcoeff);

    fockact = make_shared<Matrix>(*acoeff % *fockao * *acoeff);
    fockact_c = make_shared<Matrix>(*acoeff % *ofockao * *acoeff);
    fock   = make_shared<Matrix>(*coeffall % *fockao * *coeffall);
    fock_c = make_shared<Matrix>(*coeffall % *ofockao * *coeffall);

    // h'eff (only 1/2 exchange in the active space)
    auto fockao_p = make_shared<Fock<1>>(geom_, ofockao, ofockao->clone(), make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0))), /*store*/false, /*rhf*/false);
    fockact_p = make_shared<Matrix>(*acoeff % *fockao_p * *acoeff);
    fock_p = make_shared<Matrix>(*coeffall % *fockao_p * *coeffall);

    // h''eff (treat active orbitals as closed)
    auto fockao_h = make_shared<Fock<1>>(geom_, ofockao, nullptr, acoeff, /*store*/false, /*rhf*/true);
    fockact_h = make_shared<Matrix>(*acoeff % *fockao_h * *acoeff);
    fock_h = make_shared<Matrix>(*coeffall % *fockao_h * *coeffall);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // make K and K' matrices
  shared_ptr<Matrix> kmat;
  shared_ptr<Matrix> kmatp;
  {
    kmat = make_shared<Matrix>(*fockact_c * *rdm1);
    *kmat += Qvec(nact, nact, acoeff, /*nclosed*/0, casscf_->fci(), casscf_->fci()->rdm2(istate));

    kmatp = make_shared<Matrix>(*kmat * (-1.0));
    *kmatp += *fockact * 2.0;
  }
  shared_ptr<const Matrix> kmat2;
  shared_ptr<const Matrix> kmatp2;
  shared_ptr<const Matrix> ints2;
  {
    shared_ptr<const DFFullDist> full = casscf_->fci()->jop()->mo2e_1ext()->compute_second_transform(acoeff)->apply_J();
    // integrals (ij|kl) and <ik|jl>
    shared_ptr<const Matrix> ints = full->form_4index(full, 1.0);
    auto tmp = make_shared<Matrix>(nact*nact, nact*nact);
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(ints->data(), tmp->data(), nact, nact, nact, nact);
    ints2 = tmp;

    auto compute_kmat = [](const int nact, shared_ptr<const Matrix> rdm2, shared_ptr<const Matrix> rdm3, shared_ptr<const Matrix> fock,
                           shared_ptr<const Matrix> ints, const double sign) {
      auto out = rdm2->clone();
      // temp area
      shared_ptr<Matrix> four  = rdm2->clone();
      shared_ptr<Matrix> four2 = rdm2->clone();
      shared_ptr<Matrix> six   = rdm3->clone();

      // + (h)rdm2(i,j,k,m) * fockact_h(m,l)
      dgemm_("N", "N", nact*nact*nact, nact, nact, 1.0, rdm2->data(), nact*nact*nact, fock->data(), nact, 1.0, out->data(), nact*nact*nact);
      // + (h)rdm2(i,j,m,k) * fockact_h(m,l)
      SMITH::sort_indices<0,1,3,2,0,1,1,1>(rdm2->data(), four->data(), nact, nact, nact, nact);
      dgemm_("N", "N", nact*nact*nact, nact, nact, 1.0, four->data(), nact*nact*nact, fock->data(), nact, 0.0, four2->data(), nact*nact*nact);
      SMITH::sort_indices<0,1,3,2,1,1,1,1>(four2->data(), out->data(), nact, nact, nact, nact);
      // +/- (h)rdm2(i,j,m,n) * <mn|kl>
      dgemm_("N", "N", nact*nact, nact*nact, nact*nact, sign, rdm2->data(), nact*nact, ints->data(), nact*nact, 1.0, out->data(), nact*nact);
      // +/- (h)rdm3(i,j,x,y,l,z) * (<kx|yz>*)
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdm3->data(), six->data(), nact*nact, nact*nact, nact, nact); // (i,j,l,x,y,w)
      dgemm_("N", "T", nact*nact*nact, nact, nact*nact*nact, sign, six->data(), nact*nact*nact, ints->data(), nact, 0.0, four->data(), nact*nact*nact);
      SMITH::sort_indices<0,1,3,2,1,1,1,1>(four->data(), out->data(), nact, nact, nact, nact);
      // +/- (h)rdm3(i,j,x,k,y,z) * (<lx|yz>*)
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdm3->data(), six->data(), nact*nact, nact, nact, nact*nact); // (i,j,k,x,y,w)
      dgemm_("N", "T", nact*nact*nact, nact, nact*nact*nact, sign, six->data(), nact*nact*nact, ints->data(), nact, 1.0, out->data(), nact*nact*nact);
      return out;
    };

    kmat2  = compute_kmat(nact,  rdm2,  rdm3, fockact_c, ints2,  1.0);
    kmatp2 = compute_kmat(nact, hrdm2, hrdm3, fockact_h, ints2, -1.0);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // A matrices
  shared_ptr<Matrix> amat2 = rdm2->clone();
  {
    for (int b = 0; b != nact; ++b)
      for (int a = 0; a != nact; ++a)
        for (int bp = 0; bp != nact; ++bp)
          for (int ap = 0; ap != nact; ++ap)
            for (int c = 0; c != nact; ++c) {
              amat2->element(ap+nact*bp,a+nact*b) += fockact_p->element(c,a) * ardm2->element(bp+nact*ap,c+nact*b) - fockact_p->element(c,b) * ardm2->element(bp+nact*ap,a+nact*c);
              for (int d = 0; d != nact; ++d)
                for (int e = 0; e != nact; ++e)
                  amat2->element(ap+nact*bp,a+nact*b) += 0.5 * ints2->element(c+nact*d,e+nact*a) * (ardm3->element(id3(bp,ap,c),id3(e,d,b))
                                                                                                  + ardm3->element(id3(bp,ap,d),id3(b,c,e)))
                                                       - 0.5 * ints2->element(b+nact*c,e+nact*d) * (ardm3->element(id3(bp,ap,a),id3(e,c,d))
                                                                                                  + ardm3->element(id3(bp,ap,c),id3(d,a,e)));
            }
    shared_ptr<Matrix> tmp = amat2->copy();
    SMITH::sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), amat2->data(), nact, nact, nact, nact);
  }
  shared_ptr<Matrix> amat3 = rdm3->clone();
  {
    for (int c = 0; c != nact; ++c)
      for (int b = 0; b != nact; ++b)
        for (int a = 0; a != nact; ++a)
          for (int cp = 0; cp != nact; ++cp)
            for (int bp = 0; bp != nact; ++bp)
              for (int ap = 0; ap != nact; ++ap)
                for (int d = 0; d != nact; ++d) {
                  amat3->element(id3(ap,bp,cp),id3(a,b,c)) += fockact_p->element(d,a)*ardm3->element(id3(cp,ap,bp),id3(b,d,c))
                                                            - fockact_p->element(c,d)*ardm3->element(id3(cp,ap,bp),id3(b,a,d))
                                                            - fockact_p->element(b,d)*ardm3->element(id3(cp,ap,bp),id3(d,a,c));
                  for (int e = 0; e != nact; ++e) {
                    amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2->element(id2(c,d),id2(e,a))*ardm3->element(id3(cp,ap,bp),id3(b,d,e))
                                                          - 0.5*ints2->element(id2(d,e),id2(e,a))*ardm3->element(id3(cp,ap,bp),id3(b,d,c))
                                                          - 0.5*ints2->element(id2(c,d),id2(d,e))*ardm3->element(id3(cp,ap,bp),id3(b,a,e))
                                                          + 0.5*ints2->element(id2(b,d),id2(d,e))*ardm3->element(id3(cp,ap,bp),id3(e,a,c));
                    for (int f = 0; f != nact; ++f) {
                      amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2->element(id2(d,e),id2(f,a))*ardm4->element(id4(cp,ap,bp,b),id4(d,f,e,c))
                                                                - ints2->element(id2(d,c),id2(f,e))*ardm4->element(id4(cp,ap,bp,b),id4(d,f,a,e))
                                                                - ints2->element(id2(d,b),id2(f,e))*ardm4->element(id4(cp,ap,bp,e),id4(d,f,a,c));
                    }
                  }
                }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // B matrices
  shared_ptr<Matrix> bmat2 = make_shared<Matrix>(nact*nact*nact, nact);
  {
    for (int a = 0; a != nact; ++a)
      for (int cp = 0; cp != nact; ++cp)
        for (int bp = 0; bp != nact; ++bp)
          for (int ap = 0; ap != nact; ++ap)
            for (int c = 0; c != nact; ++c) {
              bmat2->element(id3(ap,bp,cp),a) -= (fockact_p->element(a,c)*2.0-fockact_c->element(a,c)) * ardm2->element(id2(cp,ap),id2(bp,c));
              for (int e = 0; e != nact; ++e)
                for (int f = 0; f != nact; ++f)
                  bmat2->element(id3(ap,bp,cp),a) -= ints2->element(id2(a,c),id2(e,f)) * ardm3->element(id3(cp,ap,bp),id3(e,c,f));
            }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // C matrices
  shared_ptr<Matrix> cmat2 = make_shared<Matrix>(nact, nact*nact*nact);
  {
    for (int c = 0; c != nact; ++c)
      for (int a = 0; a != nact; ++a)
        for (int b = 0; b != nact; ++b)
          for (int ap = 0; ap != nact; ++ap)
            for (int d = 0; d != nact; ++d) {
              cmat2->element(ap, id3(a,b,c)) += fockact_p->element(d,a)*ardm2->element(id2(ap,b),id2(d,c)) - fockact_p->element(c,d)*ardm2->element(id2(ap,b),id2(a,d))
                                              - fockact_p->element(b,d)*ardm2->element(id2(ap,d),id2(a,c));
              for (int e = 0; e != nact; ++e) {
                for (int f = 0; f != nact; ++f)
                  cmat2->element(ap, id3(a,b,c)) += ints2->element(id2(d,e),id2(f,a))*ardm3->element(id3(ap,b,d),id3(f,e,c))
                                                  - ints2->element(id2(d,c),id2(f,e))*ardm3->element(id3(ap,b,d),id3(f,a,e))
                                                  - ints2->element(id2(d,b),id2(f,e))*ardm3->element(id3(ap,e,d),id3(f,a,c));
                cmat2->element(ap, id3(a,b,c)) += ints2->element(id2(c,d),id2(e,a))*ardm2->element(id2(ap,b),id2(d,e))
                                             -0.5*ints2->element(id2(d,e),id2(e,a))*ardm2->element(id2(ap,b),id2(d,c))
                                             -0.5*ints2->element(id2(c,d),id2(d,e))*ardm2->element(id2(ap,b),id2(a,e))
                                             +0.5*ints2->element(id2(b,d),id2(d,e))*ardm2->element(id2(ap,e),id2(a,c));
              }
            }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // D matrices
  shared_ptr<Matrix> dmat2 = rdm2->clone();
  {
    for (int b = 0; b != nact; ++b)
      for (int a = 0; a != nact; ++a)
        for (int bp = 0; bp != nact; ++bp)
          for (int ap = 0; ap != nact; ++ap)
            for (int c = 0; c != nact; ++c) {
              dmat2->element(ap+nact*bp,a+nact*b) -= fockact_p->element(b,c) * ((a == ap ? 2.0 : 0.0) * rdm1->element(bp,c) - rdm2->element(a+nact*bp,ap+nact*c));
              dmat2->element(ap+nact*bp,a+nact*b) += fockact_p->element(a,c) * ((c == ap ? 2.0 : 0.0) * rdm1->element(bp,b) - rdm2->element(c+nact*bp,ap+nact*b));
              for (int d = 0; d != nact; ++d)
                for (int e = 0; e != nact; ++e) {
                  dmat2->element(ap+nact*bp,a+nact*b) -= 0.5 * ints2->element(c+nact*b,e+nact*d)
                           * ((a == ap ? 2.0 : 0.0) * ardm2->element(c+nact*e,bp+nact*d) - ardm3->element(id3(c,e,a),id3(ap,bp,d))
                           +  (a == ap ? 2.0 : 0.0) * ardm2->element(bp+nact*d,c+nact*e) - ardm3->element(id3(a,ap,bp),id3(d,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2->element(c+nact*e,a+nact*d)  + ardm2->element(a+nact*d,c+nact*e))
                           +  (c == ap ? 1.0 : 0.0) *((a == e  ? 2.0 : 0.0) * rdm1->element(bp, d) - ardm2->element(a+nact*e,bp+nact*d))
                           - (bp == e  ? 1.0 : 0.0) *((a == ap ? 2.0 : 0.0) * rdm1->element(c,d)   - ardm2->element(a+nact*ap,c+nact*d)));
                  dmat2->element(ap+nact*bp,a+nact*b) += 0.5 * ints2->element(c+nact*d,e+nact*a)
                           * ((d == ap ? 2.0 : 0.0) * ardm2->element(c+nact*e,bp+nact*b) - ardm3->element(id3(c,e,d),id3(ap,bp,b))
                           +  (d == ap ? 2.0 : 0.0) * ardm2->element(bp+nact*b,c+nact*e) - ardm3->element(id3(d,ap,bp),id3(b,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2->element(c+nact*e,d+nact*b)  + ardm2->element(d+nact*b,c+nact*e))
                           +  (c == ap ? 1.0 : 0.0) *((d == e  ? 2.0 : 0.0) * rdm1->element(bp, b) - ardm2->element(d+nact*e,bp+nact*b))
                           - (bp == e  ? 1.0 : 0.0) *((d == ap ? 2.0 : 0.0) * rdm1->element(c,b)   - ardm2->element(d+nact*ap,c+nact*b)));
                }
            }
    shared_ptr<Matrix> tmp = dmat2->copy();
    SMITH::sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), dmat2->data(), nact, nact, nact, nact);
  }
  shared_ptr<Matrix> dmat1 = rdm1->clone();
  {
    for (int a = 0; a != nact; ++a)
      for (int ap = 0; ap != nact; ++ap)
        for (int c = 0; c != nact; ++c) {
          dmat1->element(ap,a) += -(fockact_p->element(a,c)*2.0-fockact_c->element(a,c)) * rdm1->element(c,ap);
          for (int e = 0; e != nact; ++e)
            for (int f = 0; f != nact; ++f)
              dmat1->element(ap,a) += - ints2->element(id2(a,c),id2(e,f)) * ardm2->element(id2(ap,e),id2(c,f));
        }
  }

  Timer timer;
  // compute transformed integrals
  shared_ptr<DFDistT> fullvi;
  shared_ptr<const Matrix> fullav;
  shared_ptr<const Matrix> fullai;
  // TODO probably we want to use JKFIT for this for consistency?
  shared_ptr<const Matrix> fullaa;
  size_t memory_size;

  {
    // first compute half transformed integrals
    shared_ptr<DFHalfDist> half, halfa;
    if (abasis_.empty()) {
      half = geom_->df()->compute_half_transform(ccoeff);
      halfa = geom_->df()->compute_half_transform(acoeff);
      // used later to determine the cache size
      memory_size = half->block(0)->size() * 2;
      mpi__->broadcast(&memory_size, 1, 0);
    } else {
      auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
      auto cgeom = make_shared<Geometry>(*geom_, info, false);
      half = cgeom->df()->compute_half_transform(ccoeff);
      halfa = cgeom->df()->compute_half_transform(acoeff);
      // used later to determine the cache size
      memory_size = cgeom->df()->block(0)->size();
      mpi__->broadcast(&memory_size, 1, 0);
    }

    // second transform for virtual index and rearrange data
    {
      // this is now (naux, nvirt, nclosed), distributed by nvirt*nclosed. Always naux*nvirt block is localized to one node
      shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
      fullvi = make_shared<DFDistT>(full, dist);
    }
    {
      shared_ptr<DFFullDist> full = halfa->compute_second_transform(coeffall)->apply_J();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
      auto fullax_all = make_shared<DFDistT>(full, dist);
      shared_ptr<const Matrix> fullax = fullax_all->replicate();

      fullai = fullax->slice(0, nact*nclosed);
      fullaa = fullax->slice(nact*nclosed, nact*(nclosed+nact));
      fullav = fullax->slice(nact*(nclosed+nact), nact*(nclosed+nact+nvirt));
    }
    fullvi->discard_df();
  }
  assert(fullvi->nblocks() == 1);
  const size_t naux = fullvi->naux();

  cout << "    * 3-index integral transformation done" << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // make a list of static distribution
  const int myrank = mpi__->rank();
  vector<vector<tuple<int,int,int,int>>> tasks(mpi__->size());
  // distribution of closed-closed
  {
    StaticDist ijdist(nclosed*(nclosed+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nclosed; ++i)
        for (int j = i; j < nclosed; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j, i, /*mpitags*/-1,-1));
    }
  }
  // distribution of virt-virt (cheap as both involve active indices)
  {
    StaticDist ijdist(nvirt*(nvirt+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nvirt; ++i)
        for (int j = i; j < nvirt; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j+nclosed+nact, i+nclosed+nact, /*mpitags*/-1,-1));
    }
  }
  // distribution of closed (sort of cheap)
  {
    StaticDist ijdist(nclosed, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0; i < nclosed; ++i)
        if (i >= ijdist.start(inode) && i < ijdist.start(inode) + ijdist.size(inode))
          tasks[inode].push_back(make_tuple(i, -1, /*mpitags*/-1,-1));
    }
  }
  {
    int nmax = 0;
    for (auto& i : tasks)
      if (nmax < i.size()) nmax = i.size();
    for (auto& i : tasks) {
      const int n = i.size();
      for (int j = 0; j != nmax-n; ++j) i.push_back(make_tuple(-1,-1,-1,-1));
    }
  }
  const int nloop = tasks[0].size();

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // TODO this is identical to code in mp2.cc. Isolate

  // start communication (n fetch behind) - n is determined by memory size
  // the data is stored in a map
  map<int, shared_ptr<Matrix>> cache;
  // pair of node and set of integers
  vector<set<int>> cachetable(mpi__->size());
  vector<int> sendreqs;

  auto cache_block = [&](const int nadd, const int ndrop) {
    assert(ndrop < nadd);
    if (ndrop >= 0) {
      for (int inode = 0; inode != mpi__->size(); ++inode) {
        const int id = get<0>(tasks[inode][ndrop]);
        const int jd = get<1>(tasks[inode][ndrop]);
        // if id and jd are no longer used in the cache, delete the element
        set<int> used;
        for (int i = ndrop+1; i <= nadd; ++i) {
          used.insert(get<0>(tasks[inode][i]));
          used.insert(get<1>(tasks[inode][i]));
        }
        if (id >= 0 && id < nclosed && !used.count(id)) {
          if (inode == myrank) cache.erase(id);
          cachetable[inode].erase(id);
        }
        if (jd >= 0 && jd < nclosed && !used.count(jd)) {
          if (inode == myrank) cache.erase(jd);
          cachetable[inode].erase(jd);
        }
      }
    }
    if (nadd < nloop) {
      // issue recv requests
      auto request_one_ = [&](const int i, const int rank) {
        if (i < 0 || i >= nclosed) return -1;
        cachetable[rank].insert(i);
        int tag = -1;
        if (cache.find(i) == cache.end() && myrank == rank) {
          const int origin = fullvi->locate(0, i*nvirt);
          if (origin == myrank) {
            cache[i] = fullvi->get_slice(i*nvirt, (i+1)*nvirt).front();
          } else {
            cache[i] = make_shared<Matrix>(naux, nvirt, true);
            tag = mpi__->request_recv(cache[i]->data(), cache[i]->size(), origin, myrank*nclosed+i);
          }
        }
        return tag;
      };

      // issue send requests
      auto send_one_ = [&](const int i, const int dest) {
        // see if "i" is cached at dest
        if (i < 0 || i >= nclosed || cachetable[dest].count(i) || fullvi->locate(0, i*nvirt) != myrank)
          return -1;
        return mpi__->request_send(fullvi->data() + (i*nvirt-fullvi->bstart())*naux, nvirt*naux, dest, dest*nclosed+i);
      };

      for (int inode = 0; inode != mpi__->size(); ++inode) {
        if (inode == myrank) {
          // recieve data from other processes
          get<2>(tasks[myrank][nadd]) = request_one_(get<0>(tasks[myrank][nadd]), myrank); // receive requests
          get<3>(tasks[myrank][nadd]) = request_one_(get<1>(tasks[myrank][nadd]), myrank);
        } else {
          // send data to other processes
          const int i = send_one_(get<0>(tasks[inode][nadd]), inode); // send requests
          if (i >= 0) sendreqs.push_back(i);
          request_one_(get<0>(tasks[inode][nadd]), inode); // update cachetable
          const int j = send_one_(get<1>(tasks[inode][nadd]), inode); // send requests
          if (j >= 0) sendreqs.push_back(j);
          request_one_(get<1>(tasks[inode][nadd]), inode);
        }
      }
    }
  };
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  const int ncache = min(memory_size/(nvirt*nvirt), size_t(20));
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache_block(n, -1);

  // loop over tasks
  energy_ = 0;
double __debug = 0.0;
  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache_block(n+ncache, n-1);

    const int i = get<0>(tasks[myrank][n]);
    const int j = get<1>(tasks[myrank][n]);

    if (i < 0 && j < 0) {
      continue;
    } else if (i < nclosed && j < nclosed) {
      const int ti = get<2>(tasks[myrank][n]);
      const int tj = get<3>(tasks[myrank][n]);
      if (ti >= 0) mpi__->wait(ti);
      if (tj >= 0) mpi__->wait(tj);

      shared_ptr<const Matrix> iblock = cache.at(i);
      shared_ptr<const Matrix> jblock = cache.at(j);
      const Matrix mat(*iblock % *jblock);

      // active part
      shared_ptr<const Matrix> iablock = fullai->slice(i*nact, (i+1)*nact);
      shared_ptr<const Matrix> jablock = fullai->slice(j*nact, (j+1)*nact);
      const Matrix mat_va(*iblock % *jablock);
      const Matrix mat_av(*iablock % *jblock);
      // hole density matrix
      const Matrix mat_vaR(mat_va * *hrdm1);
      const Matrix mat_avR(*hrdm1 % mat_av);
      // K' matrix
      const Matrix mat_vaKp(mat_va * *kmatp);
      const Matrix mat_avKp(*kmatp % mat_av);

      // S(2)ij,rs sector
      const Matrix mat_aa(*iablock % *jablock);
      Matrix mat_aaR(nact, nact);
      Matrix mat_aaK(nact, nact);
      dgemv_("N", nact*nact, nact*nact, 1.0,  hrdm2->data(), nact*nact, mat_aa.data(), 1, 0.0, mat_aaR.data(), 1);
      dgemv_("N", nact*nact, nact*nact, 1.0, kmatp2->data(), nact*nact, mat_aa.data(), 1, 0.0, mat_aaK.data(), 1);
      const double norm2  = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const double denom2 = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (norm2 > norm_thresh_)
        energy_ += norm2 / (-denom2/norm2 + oeig[i]+oeig[j]);

      // TODO should thread
      // S(1)ij,r sector
      double en1 = 0.0;
      for (int v = 0; v != nvirt; ++v) {
        double norm = 0.0;
        double denom = 0.0;
        for (int a = 0; a != nact; ++a) {
          const double va = mat_va(v, a);
          const double av = mat_av(a, v);
          const double vaR = mat_vaR(v, a);
          const double avR = mat_avR(a, v);
          const double vaK = mat_vaKp(v, a);
          const double avK = mat_avKp(a, v);
          norm  += (2.0*(va*vaR + av*avR) - av*vaR - va*avR);
          denom += (2.0*(va*vaK + av*avK) - av*vaK - va*avK);
        }
        if (norm > norm_thresh_)
          en1 += norm / (-denom/norm-veig[v]+oeig[i]+oeig[j]);
      }
      if (i == j) en1 *= 0.5;
      energy_ += en1;

      // S(0)ij,rs sector
      double en = 0.0;
      for (int v = 0; v != nvirt; ++v) {
        for (int u = v+1; u < nvirt; ++u) {
          const double vu = mat(v, u);
          const double uv = mat(u, v);
          en += 2.0*(uv*uv + vu*vu - uv*vu) / (-veig[v]+oeig[i]-veig[u]+oeig[j]);
        }
        const double vv = mat(v, v);
        en += vv*vv / (-veig[v]+oeig[i]-veig[v]+oeig[j]);
      }
      if (i != j) en *= 2.0;
      energy_ += en;

    } else if (i >= nclosed+nact && j >= nclosed+nact) {
      // S(-2)rs sector
      const int iv = i-nclosed-nact;
      const int jv = j-nclosed-nact;
      shared_ptr<const Matrix> iablock = fullav->slice(iv*nact, (iv+1)*nact);
      shared_ptr<const Matrix> jablock = fullav->slice(jv*nact, (jv+1)*nact);
      Matrix mat_aa(*iablock % *jablock);
      Matrix mat_aaR(nact*nact, nact*nact);
      Matrix mat_aaK(nact*nact, nact*nact);
      dgemv_("N", nact*nact, nact*nact, 1.0,  rdm2->data(), nact*nact, mat_aa.data(), 1, 0.0, mat_aaR.data(), 1);
      dgemv_("N", nact*nact, nact*nact, 1.0, kmat2->data(), nact*nact, mat_aa.data(), 1, 0.0, mat_aaK.data(), 1);
      const double norm  = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const double denom = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (norm > norm_thresh_)
        energy_ += norm / (denom/norm - veig[iv] - veig[jv]);

    } else if (i < nclosed && j < 0) {
      // (g|vi) with i fixed
      shared_ptr<const Matrix> iblock = cache.at(i);
      // (g|ai) with i fixed
      shared_ptr<const Matrix> iablock = fullai->slice(i*nact, (i+1)*nact);
      // reordered srdm
      Matrix srdm2p(nact*nact, nact*nact);
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(srdm2->data(), srdm2p.data(), nact, nact, nact, nact);

      for (int r = 0; r != nvirt; ++r) {
        shared_ptr<const Matrix> ibr = iblock->slice(r, r+1);
        shared_ptr<const Matrix> rblock = fullav->slice(r*nact, (r+1)*nact);

        // S(-1)rs sector
        for (int s = r; s != nvirt; ++s) {
          shared_ptr<const Matrix> ibs = iblock->slice(s, s+1);
          shared_ptr<const Matrix> sblock = fullav->slice(s*nact, (s+1)*nact);
          const Matrix mat1(*ibs % *rblock); // (vi|ar) (i, r fixed)
          const Matrix mat2(*ibr % *sblock); // (vi|as) (i, s fixed)
          const Matrix mat1R(*ibs % *rblock * *rdm1); // (vi|ar) (i, r fixed)
          const Matrix mat2R(*ibr % *sblock * *rdm1); // (vi|as) (i, s fixed)
          const Matrix mat1K(*ibs % *rblock * *kmat); // (vi|ar) (i, r fixed)
          const Matrix mat2K(*ibr % *sblock * *kmat); // (vi|as) (i, s fixed)
          const double norm  = (r == s ? 1.0 : 2.0) * (mat2R.dot_product(mat2) + mat1R.dot_product(mat1) - mat2R.dot_product(mat1));
          const double denom = (r == s ? 1.0 : 2.0) * (mat2K.dot_product(mat2) + mat1K.dot_product(mat1) - mat2K.dot_product(mat1));
          if (norm > norm_thresh_)
            energy_ += norm / (denom/norm + oeig[i] - veig[r] - veig[s]);
        }

        // S(0)ir sector
        const Matrix mat1(*ibr % *fullaa); // (ir|ab)  as (1,nact*nact)
        const Matrix mat2(*rblock % *iablock); // (ra|bi) as (nact, nact)
        const Matrix mat1S(mat1 * *srdm2);
        const Matrix mat1A(mat1 * *amat2);
        const Matrix mat1Ssym(mat1S + (mat1 ^ *srdm2));
        const Matrix mat1Asym(mat1A + (mat1 ^ *amat2));
              Matrix mat2Sp(nact, nact);
              Matrix mat2D (nact, nact);
        dgemv_("N", nact*nact, nact*nact, 1.0, srdm2p.data(), nact*nact, mat2.data(), 1, 0.0, mat2Sp.data(), 1);
        dgemv_("N", nact*nact, nact*nact, 1.0, dmat2->data(), nact*nact, mat2.data(), 1, 0.0,  mat2D.data(), 1);
        const int ir = r + nclosed + nact;
        const double norm = - 2.0*mat1S.dot_product(mat1) + blas::dot_product(mat1Ssym.data(), mat1Ssym.size(), mat2.data()) + mat2Sp.dot_product(mat2)
                          + 2.0*(fock->element(ir,i) - fock_c->element(ir,i))*(fock_c->element(ir,i) + fock_h->element(ir,i))
                          + 2.0*fock_c->element(ir,i)*fock_c->element(ir,i);
        const double denom = 2.0*mat1A.dot_product(mat1) - blas::dot_product(mat1Asym.data(), mat1Asym.size(), mat2.data()) + mat2D.dot_product(mat2);
        if (norm > norm_thresh_)
          energy_ += norm / (-denom/norm + oeig[i] - veig[r]);
      }
    }
  }
mpi__->allreduce(&__debug, 1);
cout << setprecision(10) <<  __debug << endl;

  // just to double check that all the communition is done
  for (auto& i : sendreqs)
    mpi__->wait(i);
  // allreduce energy contributions
  mpi__->allreduce(&energy_, 1);

  cout << "    * assembly done" << endl << endl;
  cout << "      NEVPT2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy();
  cout << "      NEVPT2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

}
