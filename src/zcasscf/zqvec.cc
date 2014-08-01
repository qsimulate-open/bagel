//
// BAGEL - Parallel electron correlation program.
// Filename: zqvec.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/zcasscf/zqvec.h>
#include <src/rel/reloverlap.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;


ZQvec::ZQvec(const int nbasis, const int nact, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> coeff, const int nclosed,
             shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(nbasis*2, nact*2) {
  // TODO : this constructor is wrong; the 2RDM should be conjugated so the result is phase invariant. This could be corrected if the RDM class had get_conjg.

  assert(gaunt || !breit);
  if (gaunt) throw logic_error("Gaunt not implemented yet in ZQvec");
  assert((coeff->slice(nclosed*2,(nclosed+nact)*2) - *fci->jop()->coeff()).rms() < 1.0e-15);

  array<shared_ptr<const ZMatrix>,2> kcoeff = fci->kramers_coeff();
  assert(geom->nbasis()*4 == kcoeff[0]->ndim());
  assert(nact == kcoeff[0]->mdim());
  assert(nbasis*2 == coeff->mdim());


  array<list<shared_ptr<RelDFHalf>>,2> half_coulomb = fci->jop()->half_complex_coulomb();
//array<list<shared_ptr<RelDFHalf>>,2> half_gaunt   = fci->jop()->half_complex_gaunt();

  // (jr|kl) G(ji,kl) = (jr|g)(g|kl) G(ji,kl)
  //
  // [1] compute (g|kl)
  // [2] compute [g|ji] = (g|kl)*G(ji|kl)
  // [3] compute [(g|jr)->form_2index([g|ji])]*

  // in principle this is redundant, but cheap..
  // [1] compute (g|kl) on full
  unordered_map<bitset<2>, shared_ptr<const RelDFFull>> full;
  {
    array<array<shared_ptr<const Matrix>,4>,2> rocoeff;
    array<array<shared_ptr<const Matrix>,4>,2> iocoeff;
    for (int k = 0; k != 2; ++k) {
      for (int i = 0; i != 4; ++i) {
        shared_ptr<const ZMatrix> oc = kcoeff[k]->get_submatrix(i*geom->nbasis(), 0, geom->nbasis(), nact);
        rocoeff[k][i] = oc->get_real_part();
        iocoeff[k][i] = oc->get_imag_part();
      }
    }
    // full(g,kl) = J^-1_{gg'} * (g'|kl) ; complex result with kl in the krammers basis
    full = RelMOFile::compute_full(rocoeff, iocoeff, half_coulomb, /*apply_J*/false, /*apply_JJ*/true);
  }

  assert(full.size() == 4);
  unordered_map<bitset<2>, shared_ptr<RelDFFull>> full_d;
  for (auto& i : full)
    full_d.emplace(i.first, i.second->clone());

  // [2] compute [g|ji] = (g|kl)*G(ji|kl)
  // JEB : Contract 2RDM with 3idx integrals in kramers MO basis ; all indices active for 2RDM
  // JEB : shorthand -> t = target, s = source ; i.e. we are taking the "target" to be full_d, and the "source" as full
  for (auto& t : full_d) {
    for (auto& s : full) {
      bitset<4> b;
      // TODO check!
      // t^+ s^+ t s
      // JEB : read the above as if each index was a creation/annihilation operator for target and source indices
      b[3] = t.first[1]; b[2] = s.first[1]; b[1] = t.first[0]; b[0] = s.first[0];
      // JEB : take the bitset b, and return the 2rdm_av value for specified bitset
      shared_ptr<const ZRDM<2>> rdmbuf = fci->rdm2_av_kramers(b);
      // JEB : after swapping the indices order will be :
      // t^+ t s^+ s
      shared_ptr<ZRDM<2>> rdm = rdmbuf->clone();
      assert(rdm->norb() == nact);
      // JEB : the following changes the ordering G(jk|il) -> G(ji|kl) ; then contract and accumulate
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdmbuf->data(), rdm->data(), nact, nact, nact, nact);
      *t.second += *s.second->apply_2rdm(rdm);
    }
  }

  // make another transformed integral
  // compute (g|jr)
  unordered_map<bitset<1>, shared_ptr<const RelDFFull>> fullia;
  {
    array<shared_ptr<const Matrix>,4> rcoeff;
    array<shared_ptr<const Matrix>,4> icoeff;
    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> oc = coeff->get_submatrix(i*geom->nbasis(), 0, geom->nbasis(), coeff->mdim());
      rcoeff[i] = oc->get_real_part();
      icoeff[i] = oc->get_imag_part();
    }
    for (size_t t = 0; t != 2; ++t) { // JEB : loop over kramers pairs
      list<std::shared_ptr<RelDFFull>> dffull;
      for (auto& i : half_coulomb[t])
        dffull.push_back(std::make_shared<RelDFFull>(i, rcoeff, icoeff));
      DFock::factorize(dffull);
      assert(dffull.size() == 1);
      dffull.front()->scale(dffull.front()->fac()); // take care of the factor
      fullia[bitset<1>(t)] = dffull.front();
    }
  }

  // [3] compute [(g|jr)->form_2index([g|ji])]*
  unordered_map<bitset<1>, shared_ptr<ZMatrix>> qri;
  for (auto& i : fullia) { // (g|jr)
    for (auto& j : full_d) { // (g|ji)
      // JEB : match bitset values for index j
      if (i.first[0] == j.first[1]) {
        shared_ptr<ZMatrix> tmp = i.second->form_2index(j.second, 1.0, false); // JEB : contract full_d with fullia
        bitset<1> target; target[0] = j.first[0];
        if (qri.find(target) == qri.end()) {
          qri[target] = tmp;
        } else {
          *qri.at(target) += *tmp;
        }
      }
    }
  }

  // transform the active orbital to the original
  // I need overlap..
  // JEB : Transform from NaturalOrbs to standard MOs for index i
  auto overlap = make_shared<const RelOverlap>(geom);
  const ZMatView ocoeff = coeff->slice(nclosed*2, nclosed*2+nact*2);

  // JEB : conjugate needed since the above lines build up the conjugated matrix products per comment [3]
  qri[bitset<1>("0")] = qri[bitset<1>("0")]->get_conjg();
  qri[bitset<1>("1")] = qri[bitset<1>("1")]->get_conjg();

#if 0
  *this = *qri[bitset<1>("0")] * (*kcoeff[0] % *overlap * *ocoeff) + *qri[bitset<1>("1")] * (*kcoeff[1] % *overlap * *ocoeff);
#else
  this->copy_block(0,    0, ndim(), nact, qri[bitset<1>("0")]->data());
  this->copy_block(0, nact, ndim(), nact, qri[bitset<1>("1")]->data());
#endif

#if 0
  complex<double> en = 0.0;
  for (int i = 0; i != nact*2; ++i) en += element(i+nclosed*2, i) * 0.5;
  cout << setprecision(16) << " active space 2ele energy        = " << en << endl;
#endif
}


ZQvec::ZQvec(const int nbasis, const int nact, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> rcoeff, shared_ptr<const ZMatrix> acoeff, const int nclosed,
             shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(nbasis*2, nact*2) {

  assert(gaunt || !breit);
  if (gaunt) throw logic_error("Gaunt not implemented yet in ZQvec");
  assert((rcoeff->slice(nclosed*2,(nclosed+nact)*2) - *fci->jop()->coeff()).rms() < 1.0e-15);
  assert(nbasis*2 == rcoeff->mdim());

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> rrcoeff;
  array<shared_ptr<const Matrix>, 4> ircoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = acoeff->get_submatrix(i*acoeff->ndim()/4, 0, acoeff->ndim()/4, acoeff->mdim());
    shared_ptr<const ZMatrix> rc = rcoeff->get_submatrix(i*rcoeff->ndim()/4, 0, rcoeff->ndim()/4, rcoeff->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    rrcoeff[i] = rc->get_real_part();
    ircoeff[i] = rc->get_imag_part();
  }
  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom->dfs()->split_blocks();
  dfs.push_back(geom->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexr = DFock::make_half_complex(dfdists, rrcoeff, ircoeff);
  for (auto& i : half_complexr)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_factr;
  for (auto& i : half_complexr) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_factr.insert(half_complex_factr.end(), tmp.begin(), tmp.end());
  }
  half_complexr.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_factr);
  DFock::factorize(half_complex_facta);

  // (4) compute (gamma|tu)
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulltu = dffulla.front();

  // (4.5) compute (gamma|rs)
  list<shared_ptr<RelDFFull>> dffullr;
  for (auto& i : half_complex_factr)
    dffullr.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffullr);
  dffullr.front()->scale(dffullr.front()->fac()); // take care of the factor
  assert(dffullr.size() == 1);
  shared_ptr<const RelDFFull> fullrs = dffullr.front();

  // (5) form (rs|tu) where r runs fastest
  shared_ptr<const ZMatrix> rstu = fullrs->form_4index(fulltu, 1.0);
  shared_ptr<const ZMatrix> rdm2 = fci->rdm2_av()->get_conjg();
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(nbasis*2,nact*2);

  // (6) form Qrv = (rs|tu) * G(vs,tu)
  zgemm3m_("N", "N", nbasis*2, nact*2, nact*nact*nact*8, 1.0, rstu->data(), nbasis*2, rdm2->data(), nact*nact*nact*8, 0.0, out->data(), nbasis*2);

  *this = *out;
#if 0
  complex<double> en = 0.0;
  for (int i = 0; i != nact*2; ++i) en += element(i+nclosed*2, i) * 0.5;
  cout << setprecision(16) << " new active space 2ele energy        = " << en << endl;
#endif
}
