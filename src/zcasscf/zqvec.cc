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

  // in principle this is redundant, but cheap..
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
    full = RelMOFile::compute_full(rocoeff, iocoeff, half_coulomb, /*apply_J*/false, /*apply_JJ*/true);  
  }

  assert(full.size() == 4);
  unordered_map<bitset<2>, shared_ptr<RelDFFull>> full_d;
  for (auto& i : full)
    full_d.insert(make_pair(i.first, i.second->clone()));

  for (auto& t : full_d) {
    for (auto& s : full) {
      bitset<4> b;
      b[0] = s.first[0];
      b[1] = t.first[0];
      b[2] = s.first[1];
      b[3] = t.first[1];
      // t^+ s^+ t s
      shared_ptr<const ZRDM<2>> rdmbuf = fci->rdm2_av(b);
      // t^+ t s^+ s
      shared_ptr<ZRDM<2>> rdm = rdmbuf->clone();
      assert(rdm->norb() == nact);
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdmbuf->data(), rdm->data(), nact, nact, nact, nact); 
      *t.second += *s.second->apply_2rdm(rdm);
    }
  }

#if 0
  // code to check full_d
  double energy = 0.0;
  double tmp = 0.0;
  for (auto& t : full_d) {
    for (auto& s : full) {
      if (t.first == s.first) {
        array<std::shared_ptr<DFFullDist>,2> sdata = s.second->get_data();
        array<std::shared_ptr<DFFullDist>,2> tdata = t.second->get_data();
        energy += ddot_(sdata[0]->size(), sdata[0]->block(0)->get(), 1, tdata[0]->block(0)->get(), 1); 
                - ddot_(sdata[1]->size(), sdata[1]->block(0)->get(), 1, tdata[1]->block(0)->get(), 1); 
        tmp    += ddot_(sdata[0]->size(), sdata[1]->block(0)->get(), 1, tdata[0]->block(0)->get(), 1); 
                + ddot_(sdata[1]->size(), sdata[0]->block(0)->get(), 1, tdata[1]->block(0)->get(), 1); 
      }
    }
  }
  cout << setprecision(10) << energy*0.5 << endl;
  cout << setprecision(10) << tmp*0.5 << endl;
#endif

  // make another transformed integral 
  unordered_map<bitset<1>, shared_ptr<const RelDFFull>> fullia;
  {
    array<shared_ptr<const Matrix>,4> rcoeff;
    array<shared_ptr<const Matrix>,4> icoeff;
    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> oc = coeff->get_submatrix(i*geom->nbasis(), 0, geom->nbasis(), coeff->mdim());
      rcoeff[i] = oc->get_real_part();
      icoeff[i] = oc->get_imag_part();
    }
    for (size_t t = 0; t != 2; ++t) {
      list<std::shared_ptr<RelDFFull>> dffull;
      for (auto& i : half_coulomb[t])
        dffull.push_back(std::make_shared<RelDFFull>(i, rcoeff, icoeff));
      DFock::factorize(dffull);
      assert(dffull.size() == 1);
      dffull.front()->scale(dffull.front()->fac()); // take care of the factor
      fullia[bitset<1>(t)] = dffull.front();
    }
  }

  unordered_map<bitset<1>, shared_ptr<ZMatrix>> qri;
  for (auto& i : fullia) { // (g|ir)
    for (auto& j : full_d) { // (g|ij) 
      if (i.first[0] == j.first[1]) {
        shared_ptr<ZMatrix> tmp = i.second->form_2index(j.second, 1.0);
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
  auto overlap = make_shared<const RelOverlap>(geom);
  shared_ptr<const ZMatrix> ocoeff = coeff->slice(nclosed*2, nclosed*2+nact*2);

  *this = *qri[bitset<1>("0")] * (*kcoeff[0] % *overlap * *ocoeff) + *qri[bitset<1>("1")] * (*kcoeff[1] % *overlap * *ocoeff); 
}
