//
// BAGEL - Parallel electron correlation program.
// Filename: cpcasscf.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/grad/cpcasscf.h>
#include <src/math/linearRM.h>
#include <src/math/bfgs.h>
#include <src/casscf/qvec.h>

using namespace std;
using namespace bagel;

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix, Dvec>> grad, const shared_ptr<const Dvec> civ, const shared_ptr<const DFHalfDist> h,
                   const shared_ptr<const DFHalfDist> h2, const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: grad_(grad), civector_(civ), half_(h), halfjj_(h2), ref_(r), geom_(r->geom()), fci_(f) {

#if 0
  cout << "   CI vectors:" << endl;
  civector_->print(-1);
#endif
}

shared_ptr<Matrix> CPCASSCF::compute_orb_denom() const {
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nvirt = ref_->nvirt();
  const int nmobasis = ref_->coeff()->mdim();
  shared_ptr<const Matrix> acoeff = ref_->coeff()->slice(nclosed, nclosed+nact);
  shared_ptr<const Matrix> coeff = ref_->coeff();

  auto denom = make_shared<Matrix>(nmobasis, nmobasis);
  {
    denom->fill(1.0e10);
    // as in Theor Chem Acc (1997) 97:88-95
    vector<double> occup = ref_->rdm1_av()->diag();

    // inactive fock in MO basis
    auto finact = make_shared<Matrix>(*coeff % *fci_->jop()->core_fock() * *coeff);

    shared_ptr<Matrix> rdm1av = make_shared<Matrix>(nact, nact);
    copy_n(ref_->rdm1_av()->data(), rdm1av->size(), rdm1av->data());
    rdm1av->sqrt();
    rdm1av->scale(1.0/sqrt(2.0));

    auto fact_ao = make_shared<Fock<1>>(geom_, fci_->jop()->core_fock()->clone(), nullptr, make_shared<Matrix>(*acoeff * *rdm1av), /*grad*/false, /*rhf*/true);
    auto f = make_shared<Matrix>(*finact + *coeff% *fact_ao * *coeff);

    auto fact = make_shared<Qvec>(nmobasis, nact, ref_->coeff(), nclosed, fci_, ref_->rdm2_av());
    for (int i = 0; i != nact; ++i)
      daxpy_(nmobasis, occup[i], finact->element_ptr(0,nclosed+i), 1, fact->data()+i*nmobasis, 1);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nvirt; ++j)
        denom->element(j+nocc,i+nclosed) = denom->element(i+nclosed,j+nocc) = -2.0*fact->element(i,i) + 2.0*occup[i]*f->element(j+nocc, j+nocc);

    for (int i = 0; i != nclosed; ++i)
      for (int j = 0; j != nvirt; ++j)
         denom->element(j+nocc,i) = denom->element(i,j+nocc) = 4.0*f->element(j+nocc, j+nocc) - 4.0*f->element(i, i);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j)
         denom->element(j,i+nclosed) = denom->element(i+nclosed,j)
                                     = (f->element(nclosed+i,nclosed+i)*4.0-2.0*fact->element(i+nclosed,i)) - f->element(j, j)*(4.0 - 2.0*occup[i]);
  }
  return denom;
}


tuple<shared_ptr<const Matrix>, shared_ptr<const Dvec>, shared_ptr<const Matrix>>
  CPCASSCF::solve(const double zthresh, const int zmaxiter) {

  // RI determinant space
  auto detex = make_shared<Determinants>(fci_->norb(), fci_->nelea(), fci_->neleb(), false, /*mute=*/true);
  assert(fci_->norb() == ref_->nact());

  const size_t nocca = ref_->nocc();
  assert(ref_->nact() + ref_->nclosed() == nocca);

  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(0, nocca);

  // some DF vectors
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(ocoeff)->apply_J();
  shared_ptr<const DFFullDist> fullb = half->compute_second_transform(ocoeff);

  // making denominator...
  shared_ptr<PairFile<Matrix, Dvec>> denom;
  const double core_energy = geom_->nuclear_repulsion() + fci_->core_energy();
  {
    shared_ptr<Matrix> d0 = compute_orb_denom();
    for (auto& i : *d0)
      if (fabs(i) < 1.0e-10) i = 1.0;
    auto d1_tmp = make_shared<const Civec>(*fci_->denom());
    auto d1 = make_shared<Dvec>(d1_tmp, ref_->nstate());
    for (int i = 0; i != ref_->nstate(); ++i)
      *d1->data(i) -= fci_->energy(i) - core_energy;
    denom = make_shared<PairFile<Matrix, Dvec>>(d0, d1);
  }

  // BFGS update of the denominator above
  auto bfgs = make_shared<BFGS<PairFile<Matrix, Dvec>>>(denom, true);

  // CI vector
  auto source = make_shared<PairFile<Matrix, Dvec>>(*grad_);
  // antisymmetrize
  source->first()->antisymmetrize();
  source->first()->purify_redrotation(ref_->nclosed(), ref_->nact(), ref_->nvirt());

  // divide by weight
  for (int ij = 0; ij != source->second()->ij(); ++ij) {
    source->second()->data(ij)->scale(1.0/fci_->weight(ij));
  }
  // project out Civector from the gradient
  source->second()->project_out(civector_);
  auto solver = make_shared<LinearRM<PairFile<Matrix, Dvec>>>(zmaxiter, source);

  // initial guess
  shared_ptr<PairFile<Matrix, Dvec>> z = source->clone();
  z->zero();

  z = bfgs->extrapolate(source, z);
  // not needed as z->xecond() is zero
  //z->second()->project_out(civector_);

  // inverse matrix of C
  auto ovl = make_shared<const Overlap>(geom_);
  auto cinv = make_shared<const Matrix>(*ref_->coeff() % *ovl);

  cout << "  === CASSCF Z-vector iteration ===" << endl << endl;

  Timer timer;
  for (int iter = 0; iter != zmaxiter; ++iter) {
    // given z, computes sigma (before anti-symmetrization)
    shared_ptr<PairFile<Matrix, Dvec>> sigma = form_sigma(z, half, fullb, detex, cinv);
    sigma->first()->antisymmetrize();
    sigma->first()->purify_redrotation(ref_->nclosed(), ref_->nact(), ref_->nvirt());

    z = solver->compute_residual(z, sigma);

    z = bfgs->extrapolate(z, solver->civec());
    z->second()->project_out(civector_);

    cout << setw(7) <<  iter << " " << setw(20) << setprecision(14) << z->rms() << setw(15) << setprecision(2) << timer.tick() << endl;
    if (z->rms() < zthresh) break;
  }

  shared_ptr<PairFile<Matrix, Dvec>> result = solver->civec();
  shared_ptr<Matrix> xmat = form_sigma_sym(result, half, fullb, detex, cinv);

  *xmat += *grad_->first();
  xmat->symmetrize();
  xmat->scale(0.5); // due to convention
  return make_tuple(result->first(), result->second(), xmat);
}


shared_ptr<PairFile<Matrix,Dvec>>
  CPCASSCF::form_sigma(shared_ptr<const PairFile<Matrix,Dvec>> z, shared_ptr<const DFHalfDist> half,
                       shared_ptr<const DFFullDist> fullb, shared_ptr<const Determinants> detex, shared_ptr<const Matrix> cinv) const {
  const size_t nmobasis = ref_->coeff()->mdim();
  const size_t nocca = ref_->nocc();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  shared_ptr<const Matrix> z0 = z->first();
  shared_ptr<const Dvec>   z1 = z->second();
  shared_ptr<const Matrix> coeff = ref_->coeff();
  shared_ptr<RDM<1>> rdm1_av = ref_->rdm1_av()->copy();
  shared_ptr<RDM<2>> rdm2_av = ref_->rdm2_av()->copy();

  // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
  // only here we need to have det_ instead of detex
  shared_ptr<Matrix> sigmaorb = compute_amat(z1, civector_, detex);

  // computation of Atilde.
  auto cz0 = make_shared<Matrix>(*coeff * *z0);
  auto ocz0 = cz0->slice(0, nocca);
  auto cz0cinv = make_shared<Matrix>(*coeff * *z0 * *cinv);
  assert((*cz0 - *cz0cinv * *coeff).rms() < 1.0e-10);

  // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
  // (D|jx) -> (D|jS)
  {
    shared_ptr<DFFullDist> tmp0 = half->compute_second_transform(cz0cinv);
    shared_ptr<const DFHalfDist> tmp1 = geom_->df()->compute_half_transform(ocz0)->apply_J();
    tmp0->ax_plus_y(1.0, tmp1);
    shared_ptr<const DFFullDist> fulld = fullb->apply_2rdm(rdm2_av->data(), rdm1_av->data(), nclosed, nact);
    shared_ptr<const Matrix> buf = tmp0->form_2index(fulld, 2.0); // Factor of 2
    const Matrix cbuf(*coeff % *buf);
    sigmaorb->add_block(1.0, 0, 0, nmobasis, nocca, cbuf);
  }
  // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
  shared_ptr<DFFullDist> fullz = half->compute_second_transform(ocz0);
  fullz->symmetrize();
  {
    shared_ptr<const DFFullDist> tmp = fullz->apply_2rdm(rdm2_av->data(), rdm1_av->data(), nclosed, nact);
    shared_ptr<const Matrix> buf = half->form_2index(tmp, 2.0); // Factor of 2
    // mo transformation of s
    const Matrix cbuf(*coeff % *buf);
    sigmaorb->add_block(1.0, 0, 0, nmobasis, nocca, cbuf);
  }

  // one electron part...
  auto htilde = make_shared<Matrix>(*coeff % *ref_->hcore() * *cz0);
  htilde->symmetrize();
  htilde->scale(2.0);
  const Matrix cbuf(*htilde->slice(0, nocca) * *ref_->rdm1_mat());
  sigmaorb->add_block(2.0, 0, 0, nmobasis, nocca, cbuf);

  // At this point
  // htilde = Z^daggerh + hZ
  // fullb  = (D|ij)
  // fullz  = (D|ir)Z_rj + (D|rj)Z_ri

  // internal core fock operator...
  // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

  // first form 4 index
  shared_ptr<Matrix> buf = fullz->form_4index(fullb, 1.0);
  // TODO Awful code. To be updated. making the code that works in the quickest possible way
  // index swap
  unique_ptr<double[]> buf2(new double[nocca*nocca*nocca*nocca]);

  // bra ket symmetrization
  for (int i = 0; i != nocca*nocca; ++i)
    for (int j = 0; j != nocca*nocca; ++j)
      buf2[j+nocca*nocca*i] = buf->element(j, i) + buf->element(i, j);

  auto Htilde2 = make_shared<Matrix>(nact*nact, nact*nact);
  for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii)
    for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj)
      for (int k = nclosed, kk = 0; k != nocca; ++k, ++kk)
        for (int l = nclosed, ll = 0; l != nocca; ++l, ++ll)
          Htilde2->element(ll+nact*kk, jj+nact*ii) = buf2[l+nocca*(k+nocca*(j+nocca*i))];

  auto Htilde1 = make_shared<Matrix>(nact,nact, true);
  for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii) {
    for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj) {
      (*Htilde1)(jj, ii) = htilde->element(j,i);
      for (int k = 0; k != nclosed; ++k)
        (*Htilde1)(jj, ii) += 2.0*buf2[k+nocca*(k+nocca*(j+nocca*i))] - buf2[k+nocca*(i+nocca*(j+nocca*k))];
    }
  }
  // factor of 2 in the equation
  *Htilde1 *= 2.0;
  *Htilde2 *= 2.0;

  auto top = make_shared<Htilde>(ref_, 0, nact, Htilde1, Htilde2);
  vector<int> tmp(z1->ij(), 0);

  shared_ptr<Dvec> sigmaci = fci_->form_sigma(civector_, top, tmp);
  *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);

  const double core_energy = geom_->nuclear_repulsion() + fci_->core_energy();
  for (int i = 0; i != z1->ij(); ++i)
    for (int j = 0; j != z1->data(i)->size(); ++j)
      sigmaci->data(i)->data(j) -= (fci_->energy(i) - core_energy) * z1->data(i)->data(j);

  sigmaci->project_out(civector_);

  return make_shared<PairFile<Matrix, Dvec>>(sigmaorb, sigmaci);
}



// computes A matrix (scaled by 2 here)
shared_ptr<Matrix> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec, shared_ptr<const Dvec> dvec, shared_ptr<const Determinants> o) const {

  const size_t nmobasis = ref_->coeff()->mdim();

  auto amat = make_shared<Matrix>(nmobasis, nmobasis);
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  shared_ptr<const Matrix> coeff = ref_->coeff();
  shared_ptr<const Matrix> acoeff = coeff->slice(nclosed, nclosed+nact);

  // compute RDMs
  shared_ptr<const RDM<1>> rdm1t;
  shared_ptr<const RDM<2>> rdm2t;
  tie(rdm1t, rdm2t) = fci_->compute_rdm12_av_from_dvec(zvec, dvec, o);

  // symmetrize
  shared_ptr<RDM<1>> rdm1 = rdm1t->clone();
  shared_ptr<RDM<2>> rdm2 = rdm2t->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      rdm1->element(j,i) = 0.5*(rdm1t->element(j,i)+rdm1t->element(i,j));
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          rdm2->element(l,k,j,i) = 0.5*(rdm2t->element(l,k,j,i)+rdm2t->element(k,l,j,i));
  auto rdm1mat = make_shared<Matrix>(nact, nact);
  copy_n(rdm1->data(), rdm1->size(), rdm1mat->data());

  // prefactor
  const double prefactor = 2.0;

  // core Fock operator
  shared_ptr<const Matrix> core_fock = fci_->jop()->core_fock();
  amat->add_block(prefactor, 0, nclosed, nmobasis, nact, *coeff % (*core_fock * *acoeff * *rdm1mat));

  // Half transformed DF vector
  shared_ptr<const DFHalfDist> half = fci_->jop()->mo2e_1ext();
  shared_ptr<const DFFullDist> full = half->compute_second_transform(acoeff)->apply_JJ();
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(rdm2->data());
  amat->add_block(prefactor, 0, nclosed, nmobasis, nact, *coeff % *half->form_2index(fulld, 1.0));

  // additing f^z_ri contribution
  if (nclosed) {
    shared_ptr<const Matrix> aden = make_shared<Matrix>(*acoeff * *rdm1mat ^ *acoeff);
    shared_ptr<const Matrix> adenj = make_shared<Matrix>(*rdm1mat ^ *acoeff);

    shared_ptr<const Matrix> ocoeff = coeff->slice(0, nclosed);
    // coulomb
    Matrix fockz(*geom_->df()->compute_Jop(half, adenj, /*only once*/false) * *ocoeff);
    // exchange
    shared_ptr<DFFullDist> halfd = half->compute_second_transform(ocoeff);
    halfd->rotate_occ1(rdm1mat);
    fockz += *half->form_2index(halfd->apply_JJ(), -0.5);
    // add to amat
    amat->add_block(4.0, 0, 0, nmobasis, nclosed, *coeff % fockz);
  }

  return amat;
}


// as implemented in Molpro
// form sigma for calculating X (makes symmetric part correct)
shared_ptr<Matrix> CPCASSCF::form_sigma_sym(shared_ptr<const PairFile<Matrix,Dvec>> z, shared_ptr<const DFHalfDist> half,
                                            shared_ptr<const DFFullDist> fullb, shared_ptr<const Determinants> detex, shared_ptr<const Matrix> cinv) const {
  const size_t nmobasis = ref_->coeff()->mdim();
  const size_t nocca = ref_->nocc();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  shared_ptr<const Matrix> z0 = z->first();
  shared_ptr<const Dvec>   z1 = z->second();
  shared_ptr<const Matrix> coeff = ref_->coeff();
  shared_ptr<const Matrix> ccoeff = nclosed ? ref_->coeff()->slice(0, nclosed) : nullptr;
  shared_ptr<const Matrix> acoeff = ref_->coeff()->slice(nclosed, nclosed+nact);
  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(      0, nclosed+nact);

  shared_ptr<const Matrix> az0 = z0->slice(nclosed, nocca);
  shared_ptr<const Matrix> oz0 = z0->slice(0, nocca);
  shared_ptr<const Matrix> cz0  = nclosed ? z0->slice(0, nclosed) : nullptr;
  shared_ptr<const Matrix> ccz0 = nclosed ? make_shared<Matrix>(*coeff * *cz0) : nullptr;

  shared_ptr<RDM<1>> rdm1_av = ref_->rdm1_av()->copy();
  shared_ptr<RDM<2>> rdm2_av = ref_->rdm2_av()->copy();

  auto rdm1av_mat = make_shared<Matrix>(nact, nact);
  copy_n(rdm1_av->data(), rdm1_av->size(), rdm1av_mat->data());

  // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
  // only here we need to have det_ instead of detex
  shared_ptr<Matrix> sigmaorb = compute_amat(z1, civector_, detex);

  const Matrix fockinact(*coeff % *fci_->jop()->core_fock() * *coeff);

  // contributions from inactive Fock to active
  sigmaorb->add_block(2.0, 0, nclosed, nmobasis, nact, fockinact * *az0 * *rdm1av_mat);

  // TODO active Fock should be stored in somewhere
  auto weight = rdm1av_mat->copy();
  weight->sqrt();
  weight->scale(1.0/sqrt(2.0));
  const Fock<1> fact_ao(geom_, fci_->jop()->core_fock()->clone(), nullptr, make_shared<Matrix>(*acoeff * *weight), /*grad*/false, /*rhf*/true);
  const Matrix fockact(*coeff % fact_ao * *coeff);
  if (nclosed)
    sigmaorb->add_block(4.0, 0, 0, nmobasis, nclosed, (fockinact + fockact) * *cz0);

  Matrix qone(nocca, nocca);
  if (nclosed)
    qone.add_block(2.0, 0, 0, nocca, nclosed, (fockinact + fockact).get_submatrix(0, 0, nocca, nclosed));

  // TODO qvec should be stored in somewhere
  const Qvec qvec(nmobasis, nact, ref_->coeff(), nclosed, fci_, rdm2_av);
  qone.add_block(1.0, 0, nclosed, nocca, nact, qvec.get_submatrix(0, 0, nocca, nact));
  qone.add_block(1.0, 0, nclosed, nocca, nact, (*fockinact.get_submatrix(0,nclosed,nocca,nact) * *rdm1av_mat));
  qone.symmetrize();
  sigmaorb->add_block(2.0, 0, 0, nmobasis, nocca, *oz0 * qone);

  // two electron part. closed-closed
  // TODO halfc is redundant
  if (nclosed) {
    shared_ptr<DFHalfDist> halfc = geom_->df()->compute_half_transform(ccoeff);
    shared_ptr<DFHalfDist> zhalf = geom_->df()->compute_half_transform(ccz0);

    {
      shared_ptr<DFFullDist> zfull = zhalf->compute_second_transform(ccoeff)->apply_JJ();
      shared_ptr<DFFullDist> fullcj = halfc->compute_second_transform(ccoeff)->apply_JJ();
      Matrix kmat(*zhalf->form_2index(fullcj, 1.0) + *halfc->form_2index(zfull, 1.0));

      Matrix jop(*geom_->df()->compute_Jop(halfc, ccz0->transpose()) * *ccoeff * 4.0);
      sigmaorb->add_block(4.0, 0, 0, nmobasis, nclosed, *coeff % (jop - kmat));
    }
    // closed-active
    {
      shared_ptr<DFFullDist> zfull = zhalf->compute_second_transform(acoeff)->apply_JJ();
      shared_ptr<DFFullDist> fullcj = halfc->compute_second_transform(acoeff)->apply_JJ();
      Matrix kmat(*zhalf->form_2index(fullcj, 1.0) + *halfc->form_2index(zfull, 1.0));

      Matrix jop(*geom_->df()->compute_Jop(halfc, ccz0->transpose()) * *acoeff * 4.0);
      sigmaorb->add_block(2.0, 0, nclosed, nmobasis, nact, *coeff % (jop - kmat) * *rdm1av_mat);
    }
  }
  {
    shared_ptr<DFHalfDist> halfa = geom_->df()->compute_half_transform(acoeff);
    // closed-active
    if (nclosed) {
      shared_ptr<const Matrix> acz0 = make_shared<Matrix>(*coeff * *az0 * *rdm1av_mat);
      shared_ptr<DFHalfDist> zhalf = geom_->df()->compute_half_transform(acz0);
      shared_ptr<DFFullDist> zfull = zhalf->compute_second_transform(ccoeff)->apply_JJ();
      shared_ptr<DFFullDist> fullaj = halfa->compute_second_transform(ccoeff)->apply_JJ();
      Matrix kmat(*zhalf->form_2index(fullaj, 1.0) + *halfa->form_2index(zfull, 1.0));

      Matrix jop(*geom_->df()->compute_Jop(halfa, acz0->transpose()) * *ccoeff * 4.0);
      sigmaorb->add_block(2.0, 0, 0, nmobasis, nclosed, *coeff % (jop - kmat));
    }
    // active-active
    // coulomb term
    shared_ptr<const Matrix> acz0 = make_shared<Matrix>(*coeff * *az0);
    shared_ptr<DFHalfDist> zhalf = geom_->df()->compute_half_transform(acz0);
    shared_ptr<DFFullDist> fulla = halfa->compute_second_transform(acoeff);
    {
      shared_ptr<DFFullDist> fullaj = fulla->apply_JJ()->apply_2rdm(rdm2_av->data());
      shared_ptr<Matrix> jop = zhalf->form_2index(fullaj, 1.0);
      sigmaorb->add_block(2.0, 0, nclosed, nmobasis, nact, *coeff % *jop);
    }
    // exchange term
    {
      shared_ptr<RDM<2>> rdm2_av_exch = rdm2_av->clone();
      shared_ptr<DFFullDist> zfull = zhalf->compute_second_transform(acoeff)->apply_JJ();
      zfull->symmetrize();
      zfull = zfull->apply_2rdm(rdm2_av->data());
      shared_ptr<Matrix> kop = halfa->form_2index(zfull, 1.0);
      sigmaorb->add_block(2.0, 0, nclosed, nmobasis, nact, *coeff % *kop);
    }
  }

  return sigmaorb;
}
