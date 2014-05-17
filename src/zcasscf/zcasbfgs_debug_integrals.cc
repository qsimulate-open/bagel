//
// BAGEL - Parallel electron correlation program.
// Filename: zcasbfgs_debug_integrals.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson E Bates <jefferson.bates@northwestern.edu>
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

#include <src/zcasscf/zcasbfgs.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

shared_ptr<ZMatrix> ZCASBFGS::___debug___active_fock(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> rdm1) const {
  // alternative implementation of active fock matrix : returns M(a) = [ (t u|a a) - (t a|a u) ] D(t u)
  // for now we implement in the worst possible way...
  shared_ptr<ZMatrix> coefft = coeff_->slice(nclosed_*2, nocc_*2);
  assert(rdm1->mdim() == nact_*2);
 
  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> rtcoeff;
  array<shared_ptr<const Matrix>, 4> itcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coefft->get_submatrix(i*coefft->ndim()/4, 0, coefft->ndim()/4, coefft->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    rtcoeff[i] = ic->get_real_part();
    itcoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complext = DFock::make_half_complex(dfdists, rtcoeff, itcoeff);
  for (auto& i : half_complext)
    i = i->apply_J()->apply_J();

   // (3) split and factorize
   list<shared_ptr<RelDFHalf>> half_complex_factt;
   for (auto& i : half_complext) {
     list<shared_ptr<RelDFHalf>> tmp = i->split(false);
     half_complex_factt.insert(half_complex_factt.end(), tmp.begin(), tmp.end());
   }
   half_complext.clear();
 
   list<shared_ptr<RelDFHalf>> half_complex_facta;
   for (auto& i : half_complexa) {
     list<shared_ptr<RelDFHalf>> tmp = i->split(false);
     half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
   }
   half_complexa.clear();
   DFock::factorize(half_complex_factt);
   DFock::factorize(half_complex_facta);
 
   // (4) compute (gamma|xy)
   list<shared_ptr<RelDFFull>> dffulli;
   for (auto& i : half_complex_factt)
     dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
   DFock::factorize(dffulli);
   dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
   assert(dffulli.size() == 1);
   shared_ptr<const RelDFFull> fullta = dffulli.front();

   list<shared_ptr<RelDFFull>> dffulla;
   for (auto& i : half_complex_facta)
     dffulla.push_back(make_shared<RelDFFull>(i, rtcoeff, itcoeff));
   DFock::factorize(dffulla);
   dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
   assert(dffulla.size() == 1);
   shared_ptr<const RelDFFull> fullat = dffulla.front();

   // (4.5) compute (gamma|xx)
   list<shared_ptr<RelDFFull>> dffulli2;
   for (auto& i : half_complex_factt)
     dffulli2.push_back(make_shared<RelDFFull>(i, rtcoeff, itcoeff));
   DFock::factorize(dffulli2);
   dffulli2.front()->scale(dffulli2.front()->fac()); // take care of the factor
   assert(dffulli2.size() == 1);
   shared_ptr<const RelDFFull> fulltt = dffulli2.front();

   list<shared_ptr<RelDFFull>> dffulla2;
   for (auto& i : half_complex_facta)
     dffulla2.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
   DFock::factorize(dffulla2);
   dffulla2.front()->scale(dffulla2.front()->fac()); // take care of the factor
   assert(dffulla2.size() == 1);
   shared_ptr<const RelDFFull> fullaa = dffulla2.front();

  // (5) compute (tu|ab) and (tb|au)
  shared_ptr<ZMatrix> tuab = fulltt->form_4index(fullaa, 1.0);
  shared_ptr<ZMatrix> tbau = fullta->form_4index(fullat, 1.0);

  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(),coeffa->mdim());
  // (6) contract with 1RDM and accumulate
  for (int b = 0; b != coeffa->mdim(); ++b) {
    for (int a = 0; a != coeffa->mdim(); ++a) {
      for (int u = 0; u != coefft->mdim(); ++u) {
        for (int t = 0; t != coefft->mdim(); ++t) {
        // contribution from G(1,1) kramers
        (*out)(a,b) += ((*tuab)(t+coefft->mdim()*u, a+coeffa->mdim()*b) 
                       - (*tbau)(t+coefft->mdim()*b, a+coeffa->mdim()*u)) * (*rdm1)(t,u);
        }
      }
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___active_qvec_byhand(shared_ptr<const ZMatrix> coefft) const {
  // returns M(t,t) =  (t u|v w) G(v w,t u) where t is active
  // for now we implement in the worst possible way ...
  assert(coefft->mdim() == nact_*2);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coefft->mdim(),coefft->mdim());

  shared_ptr<const ZMatrix> ijkl = ___debug___all_integrals_coulomb_active(coefft);
  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();
  shared_ptr<ZMatrix> intermed1 = make_shared<ZMatrix>(*ijkl * *rdm2);

  for (int t = 0; t != coefft->mdim(); ++t) {
    for (int u = 0; u != coefft->mdim(); ++u) {
      (*out)(t,t) += (*intermed1)(t+coefft->mdim()*u,t+coefft->mdim()*u);
    }
  }

  return out;
}


/////////////////////////////////////
/// integral routines. They work. ///
/////////////////////////////////////

shared_ptr<ZMatrix> ZCASBFGS::___debug___offdiagonal_exchange_integrals(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns M(a,i) = (a i|a i) where a is an index of coeffa and i is an index of coeffi
  // for now we implement in the worst possible way ...
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }   
      
  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);
        
  // (2) half transform 
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  for (auto& i : half_complexa)
    i = i->apply_J();
    
  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }   
  half_complexa.clear();
  DFock::factorize(half_complex_facta);

  // (4) compute (gamma|xy) 
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla.front();

  // (5) form (ai|ai) where a runs fastest
  shared_ptr<const ZMatrix> aiai = fullai->form_4index(fullai, 1.0);
  for (int a = 0; a != coeffa->mdim(); ++a) {
    for (int i = 0; i != coeffi->mdim(); ++i) {
      // contribution from G(1,2)
      (*out)(a, i) += (*aiai)(a+coeffa->mdim()*i, a+coeffa->mdim()*i);
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_coulomb(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,i) = (aa|ii)
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulla = dffulla.front();

  // (5) form (aa|ii) where a runs fastest
  shared_ptr<const ZMatrix> aaii = fulla->form_4index(fulli, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a)
    for (int i = 0; i != coeffi->mdim(); ++i)
      (*out)(a, i) = (*aaii)(a+coeffa->mdim()*a, i+coeffi->mdim()*i);

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,i) = (ai|ia)
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla.front();

  // (5) form (ai|ia) where a runs fastest
  shared_ptr<const ZMatrix> aiia = fullai->form_4index(fullia, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a)
    for (int i = 0; i != coeffi->mdim(); ++i)
      (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, i+coeffi->mdim()*a); // <- only difference from the Coulomb version

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_coulomb_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,i) = (aa'|i'i)
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulla = dffulla.front();

  // (5) form (aa'|i'i) where a runs fastest
  shared_ptr<const ZMatrix> aaii = fulla->form_4index(fulli, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      (*out)(a, i) = (*aaii)(a+coeffa->mdim()*ap, ip+coeffi->mdim()*i);
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_exchange_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool off_diagonal, const bool diagonal) const {
  // returns Mat(a,i) = (ai|i'a')
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla.front();

  // (5) form (ai|ia) where a runs fastest
  shared_ptr<const ZMatrix> aiia = fullai->form_4index(fullia, 1.0);
  // need to include metric
  auto fullai_j = fullai->apply_J();
  shared_ptr<const ZMatrix> aiai = fullai_j->form_4index(fullai_j, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      // contribution from G(1,1)
      if (!off_diagonal)
        (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, ip+coeffi->mdim()*ap); // <- only difference from the Coulomb version

      // contribution from G(1,2)
      if (off_diagonal) {
        (*out)(a, i)  = (*aiai->get_conjg())(a+coeffa->mdim()*i, ap+coeffa->mdim()*ip); //  (a i|ka ki)
        (*out)(a, i) -= (*aiai->get_conjg())(a+coeffa->mdim()*ip, ap+coeffa->mdim()*i); // -(a ki|ka i)
      } else if (!diagonal) {
        (*out)(a, i) += (*aiai)(a+coeffa->mdim()*i, ap+coeffa->mdim()*ip);
        (*out)(a, i) -= (*aiai)(a+coeffa->mdim()*ip, ap+coeffa->mdim()*i);
      }
    }
  }

  return out;
}


// almost exactly the same code. Correct
shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_coulomb_active(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,j,i) = (aa|ji) where i and j are active.
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == 2*nact_);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulla = dffulla.front();

  // (5) form (aa|ii) where a runs fastest // <- only difference is here
  shared_ptr<const ZMatrix> aaii = fulla->form_4index(fulli, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim()*coeffi->mdim());
  for (int i = 0, ij = 0; i != coeffi->mdim(); ++i)
    for (int j = 0; j != coeffi->mdim(); ++j, ++ij)
      for (int a = 0; a != coeffa->mdim(); ++a)
        (*out)(a, ij) = (*aaii)(a+coeffa->mdim()*a, j+coeffi->mdim()*i);

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_exchange_active(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool with_kramers) const {
  // returns Mat(a,j,i) = (ai|ja), where i is an active index
  // with_kramers : Mat(a,j,i) = (a i|j ka) where a is an index of coeffa and i is active
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == 2*nact_);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla.front();

  // (5) form (ai|ja) where a runs fastest
  shared_ptr<const ZMatrix> aiia = fullai->form_4index(fullia, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim()*coeffi->mdim());
  for (int i = 0, ij = 0; i != coeffi->mdim(); ++i)
    for (int j = 0; j != coeffi->mdim(); ++j, ++ij)
      for (int a = 0; a != coeffa->mdim(); ++a) {
        const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
        if (with_kramers) {
          (*out)(a, ij) = (*aiia)(a+coeffa->mdim()*i, j+coeffi->mdim()*ap); // <- only difference from the Coulomb version
        } else {
          (*out)(a, ij) = (*aiia)(a+coeffa->mdim()*i, j+coeffi->mdim()*a); // <- only difference from the Coulomb version
        }
      }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_coulomb_active_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool closed_active) const {
  // returns Mat(a,t) = (a'a|t't) = (a'a|uv)G(uv,t't), where t is an active index and a is an index of coeffa
  // for closed_active : M(a,t) = (a a'|t t') = (a a'|v u)G(v u,t t')
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  // may be able to eliminate completely after hessian is confirmed
  assert(coeffi->mdim() == nact_*2);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulla = dffulla.front();

  // (5) form (ab|ji) where a runs fastest // <- only difference is here
  shared_ptr<const ZMatrix> abji = fulla->form_4index(fulli, 1.0);

  // (6) contract integrals with 2RDM
  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();
  shared_ptr<const ZMatrix> intermed1 = make_shared<ZMatrix>(*abji * *rdm2);

  // (7) form (a'a|t't) where a runs fastest
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int t = 0; t != coeffi->mdim(); ++t) {
    const int tp = (t < coeffi->mdim()/2) ? t+coeffi->mdim()/2 : t-coeffi->mdim()/2;
    for (int a = 0; a != coeffa->mdim(); ++a) {
      const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
      if (closed_active) {
        (*out)(a, t) = (*intermed1)(a+coeffa->mdim()*ap, t+coeffi->mdim()*tp);
      } else {
        (*out)(a, t) = (*intermed1)(ap+coeffa->mdim()*a, tp+coeffi->mdim()*t);
      }
    }
  }

 return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_integrals_exchange_active_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool closed_active) const {
  // returns Mat(a,t) = (a'w|va) * G(vw,t't) + (a'w|av) * G(tw,t'v) where a is an index of coeffa, and t is active
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  // may be possible to eliminate after Hessian has been confirmed
  assert(coeffi->mdim() == nact_*2);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla.front();

  // (5) form (aw|vb) where a runs fastest ; sort indices and contract with 2RDM
  shared_ptr<const ZMatrix> awvb = fullai->form_4index(fullia, 1.0);
  shared_ptr<ZMatrix> intermed1 = make_shared<ZMatrix>(coeffa->mdim()*coeffa->mdim(),coeffi->mdim()*coeffi->mdim());
  SMITH::sort_indices<0,3,2,1,0,1,1,1>(awvb->data(), intermed1->data(), coeffa->mdim(), coeffi->mdim(), coeffi->mdim(), coeffa->mdim());
  *intermed1 *= *(fci_->rdm2_av()); // stored as abtu

  // need to include metric for (aw|bv) ; contract with 2 RDM
  auto fullai_j = fullai->apply_J();
  shared_ptr<const ZMatrix> awbv = fullai_j->form_4index(fullai_j, 1.0);
  shared_ptr<ZMatrix> intermed2 = make_shared<ZMatrix>(coeffa->mdim()*coeffa->mdim(),coeffi->mdim()*coeffi->mdim()); // abwv
  SMITH::sort_indices<0,2,1,3,0,1,1,1>(awbv->data(), intermed2->data(), coeffa->mdim(), coeffi->mdim(), coeffa->mdim(), coeffi->mdim());
  shared_ptr<ZMatrix> rdmtmp = fci_->rdm2_av()->clone();
  SMITH::sort_indices<1,3,0,2,0,1,1,1>(fci_->rdm2_av()->data(), rdmtmp->data(), coeffi->mdim(), coeffi->mdim(), coeffi->mdim(), coeffi->mdim());
  *intermed2 *= *rdmtmp; // stored as abtu

  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      // contribution from G(1,1)
      if (!closed_active) 
        (*out)(a, i) = (*intermed1)(ap+coeffa->mdim()*a, ip+coeffi->mdim()*i);

      // contribution from G(1,2)
      if (closed_active) {
       (*out)(a, i) += (*intermed2->get_conjg())(ap+coeffa->mdim()*a, i+coeffi->mdim()*ip);
      } else { 
       (*out)(a, i) += (*intermed2)(ap+coeffa->mdim()*a, i+coeffi->mdim()*ip);
      }
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___all_integrals_coulomb_active(shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(i,j,k,l) = (ij|kl) where all indices are in the active space.
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == 2*nact_);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();
  DFock::factorize(half_complex_facti);

  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  // (5) form (ij|kl) where i runs fastest // <- only difference is here
  shared_ptr<ZMatrix> out = fulli->form_4index(fulli, 1.0);

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_2rdm_contraction_coulomb(shared_ptr<const ZMatrix> coeffa) const {
  // returns Mat(a,t) = (aa|vw)*(G(vw,tt)  where a is an index of coeffa, and t is active.
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  shared_ptr<ZMatrix> coefft = coeff_->slice(nclosed_*2, nocc_*2);
  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();

  // (1) compute (aa|vw) integrals
  shared_ptr<ZMatrix> maavw = ___debug___diagonal_integrals_coulomb_active(coeffa, coefft);
  assert(maavw->ndim() == coeffa->mdim());

  // (2) contract integrals with 2RDM
  shared_ptr<ZMatrix> intermed1 = make_shared<ZMatrix>(*maavw * *rdm2);
  
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coefft->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    for (int t = 0; t != coefft->mdim(); ++t) {
      (*out)(a, t) = (*intermed1)(a, t+coefft->mdim()*t) ;
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_2rdm_contraction_exchange(shared_ptr<const ZMatrix> coeffa, const bool with_kramers) const {
  // returns Mat(a,t) = (aw|va)*(G(vw,tt)  where a is an index of coeffa, and t is active.
  // with kramers : returns Mat(a,t) = (aw|v ka)*(G(vw,t kt)  where a is an index of coeffa, and t is active
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.

  shared_ptr<ZMatrix> coefft = coeff_->slice(nclosed_*2, nocc_*2);
  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();

  // (1) compute (aw|va) integrals
  shared_ptr<ZMatrix> mawva = ___debug___diagonal_integrals_exchange_active(coeffa, coefft, with_kramers); // <- only difference is here
  assert(mawva->ndim() == coeffa->mdim());

  // (2) contract integrals with 2RDM
  shared_ptr<ZMatrix> intermed1 = make_shared<ZMatrix>(*mawva * *rdm2);
  
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coefft->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    for (int t = 0; t != coefft->mdim(); ++t) {
      const int tp = (t < coefft->mdim()/2) ? t+coefft->mdim()/2 : t-coefft->mdim()/2;
      if (with_kramers) {
        (*out)(a, t) = (*intermed1)(a, t+coefft->mdim()*tp) ;
      } else {
        (*out)(a, t) = (*intermed1)(a, t+coefft->mdim()*t) ;
      }
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_1rdm_contraction_coulomb(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool with_kramers) const {
  // returns Mat(a,i) = (aa|iu) * {^{A}D}_{iu}  where a is an index of coeffa and i is active
  // for with_kramers Mat(a,i) = (a ka| ki u) * D(iu)
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_*2);
  shared_ptr<const ZMatrix> rdm1t = (transform_rdm1())->transpose();
  shared_ptr<const ZMatrix> cordm1 = make_shared<ZMatrix>(*coeffi * *rdm1t);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  array<shared_ptr<const Matrix>, 4> rrdmcoeff;
  array<shared_ptr<const Matrix>, 4> irdmcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    shared_ptr<const ZMatrix> rc = cordm1->get_submatrix(i*cordm1->ndim()/4, 0, cordm1->ndim()/4, cordm1->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
    rrdmcoeff[i] = rc->get_real_part();
    irdmcoeff[i] = rc->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute density matrix weighted (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, rrdmcoeff, irdmcoeff));
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fulli = dffulli.front();

  // (4.5) compute (gamma|xx)
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fulla = dffulla.front();

  // (5) form (aa|iu) D(iu) where a runs fastest
  shared_ptr<const ZMatrix> aaii = fulla->form_4index(fulli, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  if (with_kramers) {
    for (int a = 0; a != coeffa->mdim(); ++a) {
      const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
      for (int i = 0; i != coeffi->mdim(); ++i) {
        const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
        // G(1,1) contribution : (a ka| ki u) * D(iu)
        (*out)(a, i) = (*aaii)(a+coeffa->mdim()*ap, ip+coeffi->mdim()*i);
      }
    }
  } else {
    for (int a = 0; a != coeffa->mdim(); ++a)
      for (int i = 0; i != coeffi->mdim(); ++i)
        (*out)(a, i) = (*aaii)(a+coeffa->mdim()*a, i+coeffi->mdim()*i);
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___diagonal_1rdm_contraction_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi, const bool with_kramers) const {
  // returns Mat(a,i) = (au|ia) * {^{A}D}_{iu}  where a is an index of coeffa and i is active
  // for with_kramers, Mat(a,i) = (a i|u ka) * {^{A}D}_{u ki}
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_*2);
  shared_ptr<const ZMatrix> rdm1t = (transform_rdm1())->transpose();
  shared_ptr<const ZMatrix> cordm1 = make_shared<ZMatrix>(*coeffi * *rdm1t);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  array<shared_ptr<const Matrix>, 4> rrdmcoeff;
  array<shared_ptr<const Matrix>, 4> irdmcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    shared_ptr<const ZMatrix> rc = cordm1->get_submatrix(i*cordm1->ndim()/4, 0, cordm1->ndim()/4, cordm1->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
    rrdmcoeff[i] = rc->get_real_part();
    irdmcoeff[i] = rc->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  list<shared_ptr<RelDFHalf>> half_complex_facta;
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facti);
  DFock::factorize(half_complex_facta);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  // (4.5) compute density matrix weighted (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facta)
    dffulla.push_back(make_shared<RelDFFull>(i, rrdmcoeff, irdmcoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullardm = dffulla.front();

  // (5) form (au|ia) * {^{A}D}_{iu} where a runs fastest
  shared_ptr<const ZMatrix> aiia = fullardm->form_4index(fullia, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  if (with_kramers) {
    for (int a = 0; a != coeffa->mdim(); ++a) {
      const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
      for (int i = 0; i != coeffi->mdim(); ++i) {
        const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
        (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, ip+coeffi->mdim()*ap); // <- only difference from the Coulomb version
      }
    }
  } else {
    for (int a = 0; a != coeffa->mdim(); ++a) 
      for (int i = 0; i != coeffi->mdim(); ++i) 
        (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, i+coeffi->mdim()*a); // <- only difference from the Coulomb version
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___closed_active_offdiagonal_1rdm_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  /* returns Mat(a,i) =   [  (i ka|v a) -  (i a|v ka) ] * {^{A}D}_{v ki}
                        + [ (ki a|v ka) - (ki ka|v a) ] * {^{A}D}_{v i}
     where a is an index of coeffa and i is active
     for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  */
  assert(coeffi->mdim() == nact_*2);
  shared_ptr<const ZMatrix> rdm1t = (transform_rdm1())->transpose();
  shared_ptr<const ZMatrix> cordm1 = make_shared<ZMatrix>(*coeffi * *rdm1t);
  
  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  array<shared_ptr<const Matrix>, 4> rrdmcoeff;
  array<shared_ptr<const Matrix>, 4> irdmcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    shared_ptr<const ZMatrix> rc = cordm1->get_submatrix(i*cordm1->ndim()/4, 0, cordm1->ndim()/4, cordm1->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
    rrdmcoeff[i] = rc->get_real_part();
    irdmcoeff[i] = rc->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform ; traditional on 2nd quantities ; J^-1 applied to RDM weighted
  list<shared_ptr<RelDFHalf>> half_complexi  = DFock::make_half_complex(dfdists, rrdmcoeff, irdmcoeff);
  list<shared_ptr<RelDFHalf>> half_complexi2 = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti; // density matrix weighted
  list<shared_ptr<RelDFHalf>> half_complex_facti2; // traditional
  for (auto& i : half_complexi2) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti2.insert(half_complex_facti2.end(), tmp.begin(), tmp.end());
  }
  half_complexi2.clear();
  DFock::factorize(half_complex_facti2);
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();
  DFock::factorize(half_complex_facti);

  // (4) compute density matrix weighted (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facti)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullrdma = dffulla.front();

  // (4.5) compute traditional (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulla2;
  for (auto& i : half_complex_facti2)
    dffulla2.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla2);
  dffulla2.front()->scale(dffulla2.front()->fac()); // take care of the factor
  assert(dffulla2.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulla2.front();

  shared_ptr<const ZMatrix> iaia = fullia->form_4index(fullrdma, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      // G(1,2) contributions
      (*out)(a, i)  = (*iaia)(i+coeffi->mdim()*ap, ip+coeffi->mdim()*a);    //   (i ka|v a) D(v ki)
      (*out)(a, i) -= (*iaia)(i+coeffi->mdim()*a,  ip+coeffi->mdim()*ap);   // - (i a|v ka) D(v ki)
      (*out)(a, i) += (*iaia)(ip+coeffi->mdim()*a, i+coeffi->mdim()*ap);    // + (ki a|v ka) D(v i)
      (*out)(a, i) -= (*iaia)(ip+coeffi->mdim()*ap, i+coeffi->mdim()*a);    // - (ki ka|v a) D(v i)
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___closed_active_offdiagonal_2rdm_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,i) = (u ka|v a) * G(u i, v ki)   where a is an index of coeffa and i is active
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_*2);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform
  list<shared_ptr<RelDFHalf>> half_complexi = DFock::make_half_complex(dfdists, ricoeff, iicoeff);
  for (auto& i : half_complexi)
    i = i->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti;
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();

  DFock::factorize(half_complex_facti);


  // (4) compute (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulli;
  for (auto& i : half_complex_facti)
    dffulli.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff)); // <- only difference from the Coulomb version
  DFock::factorize(dffulli);
  dffulli.front()->scale(dffulli.front()->fac()); // take care of the factor
  assert(dffulli.size() == 1);
  shared_ptr<const RelDFFull> fullia = dffulli.front();

  // (5) compute (i a|j b)
  shared_ptr<ZMatrix> iajb = fullia->form_4index(fullia, 1.0);
  shared_ptr<ZMatrix> intermed1 = make_shared<ZMatrix>(coeffa->mdim()*coeffa->mdim(),nact_*nact_*4);
  SMITH::sort_indices<1,3,0,2,0,1,1,1>(iajb->data(), intermed1->data(), coeffi->mdim(), coeffa->mdim(), coeffi->mdim(), coeffa->mdim()); // sorted to abij
  shared_ptr<ZMatrix> rdm2tmp = fci_->rdm2_av()->clone();
  SMITH::sort_indices<0,2,1,3,0,1,1,1>(fci_->rdm2_av()->data(), rdm2tmp->data(), coeffi->mdim(), coeffi->mdim(), coeffi->mdim(), coeffi->mdim());

  *intermed1 *= *(rdm2tmp); 

  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      // contribution from G(1,2)
      (*out)(a, i) = (*intermed1)(a+coeffa->mdim()*ap, ip+coeffi->mdim()*i);
    }
  }

  return out;
}


shared_ptr<ZMatrix> ZCASBFGS::___debug___closed_active_diagonal_1rdm_contraction_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,i) =  (a i|u ka) * D(u ki) where a is an index of coeffa and i is active (a t|u ka) * {^{A}D}_{u kt}
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_*2);
  shared_ptr<const ZMatrix> rdm1t = (transform_rdm1())->transpose();
  shared_ptr<const ZMatrix> cordm1 = make_shared<ZMatrix>(*coeffi * *rdm1t);

  // (1) Sepeate real and imaginary parts for pcoeff
  array<shared_ptr<const Matrix>, 4> racoeff;
  array<shared_ptr<const Matrix>, 4> iacoeff;
  array<shared_ptr<const Matrix>, 4> ricoeff;
  array<shared_ptr<const Matrix>, 4> iicoeff;
  array<shared_ptr<const Matrix>, 4> rrdmcoeff;
  array<shared_ptr<const Matrix>, 4> irdmcoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> ac = coeffa->get_submatrix(i*coeffa->ndim()/4, 0, coeffa->ndim()/4, coeffa->mdim());
    shared_ptr<const ZMatrix> ic = coeffi->get_submatrix(i*coeffi->ndim()/4, 0, coeffi->ndim()/4, coeffi->mdim());
    shared_ptr<const ZMatrix> rc = cordm1->get_submatrix(i*cordm1->ndim()/4, 0, cordm1->ndim()/4, cordm1->mdim());
    racoeff[i] = ac->get_real_part();
    iacoeff[i] = ac->get_imag_part();
    ricoeff[i] = ic->get_real_part();
    iicoeff[i] = ic->get_imag_part();
    rrdmcoeff[i] = rc->get_real_part();
    irdmcoeff[i] = rc->get_imag_part();
  }

  // (1.5) dfdists
  vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
  dfs.push_back(geom_->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) half transform ; traditional on 2nd quantities ; J^-1 applied to RDM weighted
  list<shared_ptr<RelDFHalf>> half_complexi  = DFock::make_half_complex(dfdists, rrdmcoeff, irdmcoeff);
  list<shared_ptr<RelDFHalf>> half_complexa = DFock::make_half_complex(dfdists, racoeff, iacoeff);
  for (auto& i : half_complexi)
    i = i->apply_J()->apply_J();

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_facti; // density matrix weighted
  list<shared_ptr<RelDFHalf>> half_complex_facta; // traditional
  for (auto& i : half_complexa) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facta.insert(half_complex_facta.end(), tmp.begin(), tmp.end());
  }
  half_complexa.clear();
  DFock::factorize(half_complex_facta);
  for (auto& i : half_complexi) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_facti.insert(half_complex_facti.end(), tmp.begin(), tmp.end());
  }
  half_complexi.clear();
  DFock::factorize(half_complex_facti);

  // (4) compute  density matrix weighted (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulla;
  for (auto& i : half_complex_facti)
    dffulla.push_back(make_shared<RelDFFull>(i, racoeff, iacoeff));
  DFock::factorize(dffulla);
  dffulla.front()->scale(dffulla.front()->fac()); // take care of the factor
  assert(dffulla.size() == 1);
  shared_ptr<const RelDFFull> fullrdma = dffulla.front();

  // (4.5) compute traditional (gamma|xy)
  list<shared_ptr<RelDFFull>> dffulla2;
  for (auto& i : half_complex_facta)
    dffulla2.push_back(make_shared<RelDFFull>(i, ricoeff, iicoeff));
  DFock::factorize(dffulla2);
  dffulla2.front()->scale(dffulla2.front()->fac()); // take care of the factor
  assert(dffulla2.size() == 1);
  shared_ptr<const RelDFFull> fullai = dffulla2.front();

  // (5) form (ai|u ka) * {^{A}D}_{u ki} where a runs fastest
  shared_ptr<const ZMatrix> aiia = fullai->form_4index(fullrdma, 1.0);
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim());
  for (int a = 0; a != coeffa->mdim(); ++a) {
    const int ap = (a < coeffa->mdim()/2) ? a+coeffa->mdim()/2 : a-coeffa->mdim()/2;
    for (int i = 0; i != coeffi->mdim(); ++i) {
      const int ip = (i < coeffi->mdim()/2) ? i+coeffi->mdim()/2 : i-coeffi->mdim()/2;
      (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, ip+coeffi->mdim()*ap); // <- only difference from the Coulomb version
    }
  }

  return out;
}
