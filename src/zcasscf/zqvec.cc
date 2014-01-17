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

  assert(gaunt || !breit);
  if (gaunt) throw logic_error("Gaunt not implemented yet in ZQvec");

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
    array<array<shared_ptr<const Matrix>,4>,2> rocoeff; // JEB : array to hold real part of the four component AO basis with two blocks for the kramers basis
    array<array<shared_ptr<const Matrix>,4>,2> iocoeff; // JEB : array to hold imag part of the four component AO basis with two blocks for the kramers basis
    for (int k = 0; k != 2; ++k) {
      for (int i = 0; i != 4; ++i) {
        shared_ptr<const ZMatrix> oc = kcoeff[k]->get_submatrix(i*geom->nbasis(), 0, geom->nbasis(), nact);
        // JEB : double operator construct means take the kth element of rocoeff, then extract the ith element of the result.
        //       in this case then, the k-index runs over kramers pairs, and the i index runs over the 4 components of the AO basis
        //      e.g. for k=0, i=2 the result would be S^+ for the first kramers pair (in addition to the real and imaginary parts)
        rocoeff[k][i] = oc->get_real_part(); 
        iocoeff[k][i] = oc->get_imag_part();
      }
    }
    // JEB : result of compute_full is to return a three idx ERI in the kramers MO basis (with re and im parts), e.g. a transformation of Eq. (36) in KS_JCP13 to all MO
    // JEB : full(g,kl) = J^-1_{gg'} * (g'|kl) ; complex result with kl in the krammers basis
    full = RelMOFile::compute_full(rocoeff, iocoeff, half_coulomb, /*apply_J*/false, /*apply_JJ*/true); 
  } 

  // JEB : copy over the size,dimensions,bitsets of full to full_d, but no values
  assert(full.size() == 4); // JEB: full is size 4 from kramers*complex quantity
  unordered_map<bitset<2>, shared_ptr<RelDFFull>> full_d;
  for (auto& i : full)
    full_d.insert(make_pair(i.first, i.second->clone())); 

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
      shared_ptr<const ZRDM<2>> rdmbuf = fci->rdm2_av(b); 
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
        // JEB : if this is the first time the kramers index apperas, put the value of form_2indx on qri, otherwise add onto previous value
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
  shared_ptr<const ZMatrix> ocoeff = coeff->slice(nclosed*2, nclosed*2+nact*2);

  // JEB : conjugate needed since the above lines build up the conjugated matrix products per comment [3]
  qri[bitset<1>("0")] = qri[bitset<1>("0")]->get_conjg();
  qri[bitset<1>("1")] = qri[bitset<1>("1")]->get_conjg();

  *this = *qri[bitset<1>("0")] * (*kcoeff[0] % *overlap * *ocoeff) + *qri[bitset<1>("1")] * (*kcoeff[1] % *overlap * *ocoeff);

#if 0
  complex<double> en = 0.0;
  for (int i = 0; i != nact*2; ++i) en += element(i+nclosed*2, i) * 0.5;
  cout << setprecision(10) << en << endl;
#endif
}
