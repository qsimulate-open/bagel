//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_debug.cc
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

#include <src/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;

static const int morbital = 2;
static const int norbital = 3;

void ZCASSCF::___debug___orbital_rotation(const bool kramers) {
  // currently transforming n and mth orbitals;
  const double angle = idata_->get<double>("debugrot", 0.0);
  if (angle == 0.0) return;

  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << "perturbing an orbital by : " << setprecision(5) << angle << endl;
  cout << "kramers adaptation       : " << (kramers ? "on" : "off") << endl;

  shared_ptr<ZRotFile> atmp = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2, /*superci*/false);

  // currently ++ and -- blocks
  atmp->ele_vc(morbital, norbital) = angle;
  if (kramers)
    atmp->ele_vc(nvirt_+morbital, nclosed_+norbital) = angle; // = conj(angle)

  shared_ptr<ZMatrix> amattmp = atmp->unpack<ZMatrix>();

  // multiply -1 from the formula. multiply -i to make amat hermite (will be compensated)
  *amattmp *= -1.0 * complex<double>(0.0, -1.0);

  // restore the matrix from RotFile
  unique_ptr<double[]> eig(new double[amattmp->ndim()]);
  amattmp->diagonalize(eig.get());

  auto amattmp_sav = amattmp->copy();
  for (int i = 0; i != amattmp->ndim(); ++i) {
    complex<double> ex = exp(complex<double>(0.0, eig[i]));
    for_each(amattmp->element_ptr(0,i), amattmp->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
  }

  auto expatmp = make_shared<ZMatrix>(*amattmp ^ *amattmp_sav);

  if (kramers)
    kramers_adapt(expatmp);

  shared_ptr<const ZMatrix> pco = make_shared<const ZMatrix>(*coeff_ * *expatmp);
  auto tt = make_shared<ZMatrix>(*coeff_ - *pco);

  cout << "norm delta(coeff) = " << setprecision(8) << tt->norm() << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;

  coeff_ = pco;
}


void ZCASSCF::___debug___print_gradient(shared_ptr<const ZRotFile> grad, const bool with_kramers) const {
  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << "orbital gradient" << endl;

  // currently ++ and -- blocks
  const complex<double> gradient = grad->ele_vc(morbital, norbital)
                                 + (with_kramers ? grad->ele_vc(nvirt_+morbital, nclosed_+norbital) : complex<double>(0.0));

  cout << setprecision(10) << gradient << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;
}


void ZCASSCF::___debug___compute_hessian(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, const bool with_kramers) const {
  shared_ptr<const ZMatrix> coeffa = coeff_->slice(nocc_*2, coeff_->mdim());
  shared_ptr<const ZMatrix> coeffi = coeff_->slice(0, nclosed_*2);
  shared_ptr<ZMatrix> maaii = ___debug___diagonal_integrals_coulomb(coeffa, coeffi);
  shared_ptr<ZMatrix> maiia = ___debug___diagonal_integrals_exchange(coeffa, coeffi);
  *maiia -= *maaii;

  for (int i = 0; i != nclosed_*2; ++i) {
    for (int a = 0; a != nvirt_*2; ++a) {
      const int na = a + nocc_*2;
      (*maiia)(a, i) += (*cfock)(na,na) + (*afock)(na,na) - (*cfock)(i,i) - (*afock)(i,i);
    }
  }

  if (with_kramers) {
    shared_ptr<ZMatrix> kaaii = ___debug___diagonal_integrals_coulomb_kramers(coeffa, coeffi);
    shared_ptr<ZMatrix> kaiia = ___debug___diagonal_integrals_exchange_kramers(coeffa, coeffi);
    *maiia += *kaiia;
    *maiia -= *kaaii; // this appears zero for (at least) coulomb integrals due to symmetry
    *maiia += *maiia->get_conjg(); // due to ++ <-> -- and +- <-> -+ symmetry
  }

  if (nact_) {
    double fcienergy = ___debug___recompute_fci_energy(cfock->get_submatrix(nclosed_*2,nclosed_*2,nact_*2,nact_*2));
    cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
    cout << "recomputed FCI energy" << endl;
    cout << setprecision(12) << fcienergy << endl;
    cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;
  }

  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << "diagonal hessian value" << endl;
  cout << setprecision(10) << (*maiia)(morbital, norbital)*2.0 << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl << endl;
}


double ZCASSCF::___debug___recompute_fci_energy(shared_ptr<const ZMatrix> cfock) const {
  // returns FCI energy ; requires core fock matrix for active orbitals as input
  assert(cfock->ndim() == 2*nact_ && cfock->ndim() == cfock->mdim());

  shared_ptr<const ZMatrix> rdm2 = fci_->rdm2_av();
  shared_ptr<const ZMatrix> rdm1 = transform_rdm1();
  shared_ptr<ZMatrix> ijkl = ___debug___all_integrals_coulomb_active(coeff_->slice(nclosed_*2, nclosed_*2 + nact_*2));

  double twoelen = (rdm2->dot_product(ijkl)).real();
  double oneelen = (rdm1->dot_product(cfock)).real();
  double tote = fci_->core_energy()  + geom_->nuclear_repulsion() + oneelen + 0.5*twoelen;

  return tote;
}

/////////////////////////////////////
/// integral routines. They work. ///
/////////////////////////////////////

shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_coulomb(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
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


shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_exchange(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
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



shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_coulomb_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
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

  // (5) form (aa|ii) where a runs fastest
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


shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_exchange_kramers(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
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
      (*out)(a, i) = (*aiia)(a+coeffa->mdim()*i, ip+coeffi->mdim()*ap); // <- only difference from the Coulomb version

      // contribution from G(1,2)
      (*out)(a, i) += (*aiai)(a+coeffa->mdim()*i, ap+coeffa->mdim()*ip);
      (*out)(a, i) -= (*aiai)(a+coeffa->mdim()*ip, ap+coeffa->mdim()*i);
    }
  }

  return out;
}



// almost exactly the same code. Correct
shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_coulomb_active(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,j,i) = (aa|ji) where i and j are active.
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_);

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


shared_ptr<ZMatrix> ZCASSCF::___debug___diagonal_integrals_exchange_active(shared_ptr<const ZMatrix> coeffa, shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(a,j,i) = (ai|ja), where i is an active index
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_);

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
  shared_ptr<ZMatrix> out = make_shared<ZMatrix>(coeffa->mdim(), coeffi->mdim()*coeffi->mdim());
  for (int i = 0, ij = 0; i != coeffi->mdim(); ++i)
    for (int j = 0; j != coeffi->mdim(); ++j, ++ij)
      for (int a = 0; a != coeffa->mdim(); ++a)
        (*out)(a, ij) = (*aiia)(a+coeffa->mdim()*i, j+coeffi->mdim()*a); // <- only difference from the Coulomb version

  return out;
}


shared_ptr<ZMatrix> ZCASSCF::___debug___all_integrals_coulomb_active(shared_ptr<const ZMatrix> coeffi) const {
  // returns Mat(i,j,k,l) = (ij|kl) where all indices are in the active space.
  // for the time being, we implement it in the worst possible way... to be updated to make it efficient.
  assert(coeffi->mdim() == nact_);

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
