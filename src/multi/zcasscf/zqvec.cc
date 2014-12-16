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

#include <src/multi/zcasscf/zqvec.h>
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;


ZQvec::ZQvec(const int nbasis, const int nact, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> rcoeff, shared_ptr<const ZMatrix> acoeff, const int nclosed,
             shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(nbasis, nact*2) {

  assert(gaunt || !breit);
  assert((*acoeff - *fci->jop()->coeff()).rms() < 1.0e-15);
  assert(nbasis == rcoeff->mdim());

  // (1) Sepeate real and imaginary parts for coeffs
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
  shared_ptr<ZMatrix> out;
  auto compute = [&racoeff, &iacoeff, &rrcoeff, &ircoeff, &geom, &fci, &nact] (shared_ptr<ZMatrix>& out, const bool gaunt, const bool breit = false) {
   // purpose : compute (ps|tu) Gamma_{tu,ws} for coulomb, gaunt and breit ; a factor -1.0 (-0.5) and different df objects are needed for gaunt (breit)
   // (1.5) dfdists
   vector<shared_ptr<const DFDist>> dfs;
   if (!gaunt) {
     dfs = geom->dfs()->split_blocks();
     dfs.push_back(geom->df());
   } else {
     dfs = geom->dfsl()->split_blocks();
   }
   list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, gaunt);

   // (2) half transform ; half_complexa = J^{-1/2}_{gamma,delta}(delta|i mu) ; half_complexa_breit = (delta|i mu)
   list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
   list<shared_ptr<RelDFHalf>> half_complexa_breit;
   if (breit) half_complexa_breit = DFock::make_half_complex(dfdists, racoeff, iacoeff);
   for (auto& i : half_complexa)
     i = i->apply_J();

   // (3) split and factorize ; half_complex_facta = J^{-1/2}_{gamma,delta}(delta|i mu) ; half_complex_facta_breit = (delta|i mu)
   list<shared_ptr<RelDFHalf>> half_complex_facta, half_complex_facta2, half_complex_facta_breit;
   for (auto& i : half_complexa) {
     list<shared_ptr<RelDFHalf>> tmp = i->split(false);
     half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
   }
   half_complexa.clear();
   DFock::factorize(half_complex_facta);
   if (breit) {
     for (auto& i : half_complexa_breit) {
       list<shared_ptr<RelDFHalf>> tmp = i->split(false);
       half_complex_facta_breit.insert(half_complex_facta_breit.end(), tmp.begin(), tmp.end());
     }
     half_complexa_breit.clear();
     DFock::factorize(half_complex_facta_breit);
   }

   // (3.5) multiply breit 2index metric ; half_complex_facta2 = ( B^{nu nu'} J^{-1/2} )_{gamma,delta} (delta|i mu)
   if (breit) {
     // TODO Not the best implementation -- one could avoid apply_J to half-transformed objects
     auto breitint = make_shared<BreitInt>(geom);
     list<shared_ptr<Breit2Index>> breit_2index;
     for (int i = 0; i != breitint->Nblocks(); ++i) {
       breit_2index.push_back(make_shared<Breit2Index>(breitint->index(i), breitint->data(i), geom->df()->data2()));
       if (breitint->not_diagonal(i))
         breit_2index.push_back(breit_2index.back()->cross());
     }
     for (auto& i : half_complex_facta_breit)
       half_complex_facta2.push_back(i->apply_J());

     for (auto& i : half_complex_facta_breit)
       for (auto& j : breit_2index)
         if (i->alpha_matches(j)) {
           half_complex_facta2.push_back(i->apply_J()->multiply_breit2index(j));
           DFock::factorize(half_complex_facta2);
         }
   }

   // (4) compute J^{-1/2}_{delta,gamma} (gamma|tu)
   list<shared_ptr<RelDFFull>> dffulla;
   for (auto& i : half_complex_facta)
     dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
   DFock::factorize(dffulla);
   dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
   assert(dffulla.size() == 1);
   shared_ptr<const RelDFFull> fulltu = dffulla.front();

   // (4.25) compute ( B^{nu nu'} J^{-1/2} )_{delta,gamma}(gamma|tu)
   list<shared_ptr<RelDFFull>> dffulla_breit;
   shared_ptr<const RelDFFull> fulltu_breit;
   if (breit) {
     for (auto& i : half_complex_facta2)
       dffulla_breit.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
     DFock::factorize(dffulla_breit);
     dffulla_breit.front()->scale(dffulla_breit.front()->fac()); // take care of the factor
     assert(dffulla_breit.size() == 1);
     fulltu_breit = dffulla_breit.front();
   }

   // (4.5) compute J^{-1/2}_{delta,gamma} (gamma|rs)
   list<shared_ptr<RelDFFull>> dffullr;
   for (auto& i : half_complex_facta)
     dffullr.push_back(make_shared<RelDFFull>(i, rrcoeff, ircoeff)); // <- only difference from the Coulomb version
   DFock::factorize(dffullr);
   dffullr.front()->scale(dffullr.front()->fac()); // take care of the factor
   assert(dffullr.size() == 1);
   shared_ptr<const RelDFFull> fullrs = dffullr.front();

   // (5) form (rs|tu)*G(vs,tu) where r runs fastest
   auto rdm2_av = make_shared<ZRDM<2>>(nact*2);
   copy_n(fci->rdm2_av()->data(), nact*nact*nact*nact*16, rdm2_av->data());

   shared_ptr<const RelDFFull> fulltu_d = !breit ? fulltu->apply_2rdm(rdm2_av) : fulltu_breit->apply_2rdm(rdm2_av);
   const double gscale = breit ? -0.5 : -1.0;
   if (!gaunt)
     out = fullrs->form_2index(fulltu_d, 1.0, false);
   else if (!breit)
     *out += *fullrs->form_2index(fulltu_d, gscale, false);
   else {
     // symmetrization was used in ZFCI compute_mo2e, but following that procedure leads to large errors here
     // TODO : ensure that symmetrization is not needed
     *out += *fullrs->form_2index(fulltu_d, gscale, false);
   }
  };
  compute(out, false);
  if (gaunt)
    compute(out, gaunt, breit);

  *this = *out;

#if 0
  complex<double> en = 0.0;
  for (int i = 0; i != nact*2; ++i) en += element(i+nclosed*2, i) * 0.5;
  cout << setprecision(16) << " new active space 2ele energy        = " << en << endl;
#endif
}
