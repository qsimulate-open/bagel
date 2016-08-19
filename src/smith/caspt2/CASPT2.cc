//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <iostream>
#include <iomanip>
#include <src/smith/caspt2/CASPT2.h>
#include <src/util/math/linearRM.h>
#include <src/smith/caspt2/MSCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = ref->ciwfn()->nstates();

  // MS-CASPT2: t2 and s as MultiTensor (t2all, sall)
  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp)
      j = init_amplitude();
    t2all_.push_back(tmp);

    auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp2)
      j = init_residual();
    sall_.push_back(tmp2);
    rall_.push_back(tmp2->clone());
  }
  energy_.resize(nstates_);
  pt2energy_.resize(nstates_);
}


void CASPT2::CASPT2::solve() {
  Timer timer;
  print_iteration();

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  for (int istate = 0; istate != nstates_; ++istate) { //K states
    t2all_[istate]->fac(istate) = 0.0;
    sall_[istate]->fac(istate)  = 0.0;

    for (int jst=0; jst != nstates_; ++jst) { // <jst|
      if (info_->sssr() && jst != istate)
        continue;
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
      while(!sourceq->done())
        sourceq->next_compute();
    }
  }

  // solve linear equation for t amplitudes
  t2all_ = solve_linear(sall_, t2all_);
  timer.tick_print("CASPT2 energy evaluation");
  cout << endl;

  vector<double> shift_correction(nstates_*nstates_);
  for (int istate = 0; istate != nstates_; ++istate) {
    if (info_->shift() == 0.0) {
      pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate);
      cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] <<endl;
    } else {
      // will be used in normq
      n = init_residual();
      double norm = 0.0;
      for (int jstate = 0; jstate != nstates_; ++jstate) {
        if (istate != jstate && info_->shift_diag())
          continue;
        double nn = 0.0;
        for (int jst = 0; jst != nstates_; ++jst) { // bra
          for (int ist = 0; ist != nstates_; ++ist) { // ket
            if (info_->sssr() && (jst != istate || ist != istate))
              continue;
            set_rdm(jst, ist);
            t2 = t2all_[istate]->at(ist);
            shared_ptr<Queue> normq = make_normq(true, jst == ist);
            while (!normq->done())
              normq->next_compute();
            nn += dot_product_transpose(n, t2all_[istate]->at(jst));
          }
        }
        shift_correction[istate + nstates_*jstate] = nn;
        if (jstate == istate) norm = nn;
      }

      pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate) - info_->shift()*norm;
      cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
      cout << "        w/o shift correction  " << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) <<endl;
      cout <<endl;
    }
  }

  // MS-CASPT2
  if (info_->shift() && info_->do_ms())
    cout << "    MS-CASPT2:  Applying levelshift correction to " << (info_->shift_diag() ? "diagonal" : "all" ) <<  " elements of the effective Hamiltonian."  << endl << endl;
  if (info_->do_ms() && info_->sssr())
    for (int istate = 0; istate != nstates_; ++istate) //K states
      for (int jst=0; jst != nstates_; ++jst) // <jst|
        if (info_->sssr() && jst != istate) {
          set_rdm(jst, istate);
          s = sall_[istate]->at(jst);
          shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
          while(!sourceq->done())
            sourceq->next_compute();
        }

  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<Matrix>(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          (*heff_)(ist, ist) = pt2energy_[ist];
        } else if (ist < jst) {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > ]
          (*heff_)(jst, ist) = 0.5*(detail::real(dot_product_transpose(sall_[ist], t2all_[jst]))
                                  + detail::real(dot_product_transpose(sall_[jst], t2all_[ist])))
                             + (*eref_)(jst, ist);
          if (!info_->shift_diag())
            (*heff_)(jst, ist) -= shift_correction[jst+nstates_*ist]*info_->shift();
          (*heff_)(ist, jst) = (*heff_)(jst, ist);
        }
      }
    }

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASPT2 Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    heff_->diagonalize(eig);
    copy_n(eig.data(), nstates_, pt2energy_.data());

    // print out the eigen vector
    cout << endl;
    cout << "    * MS-CASPT2 rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate)
      cout << "    * MS-CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    cout << endl << endl;
  } else {
    heff_ = make_shared<Matrix>(1,1);
    heff_->element(0,0) = 1.0;
  }
  energy_ = pt2energy_;
}


// function to solve linear equation
vector<shared_ptr<MultiTensor_<double>>> CASPT2::CASPT2::solve_linear(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  Timer mtimer;
  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;
  for (int i = 0; i != nstates_; ++i) {  // K states
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i] - info_->shift();
    energy_[i] = 0.0;
    // set guess vector
    t[i]->zero();
    if (s[i]->rms() < 1.0e-15) {
      print_iteration(0, 0.0, 0.0, mtimer.tick());
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      update_amplitude(t[i], s[i]);
    }

    auto solver = make_shared<LinearRM<MultiTensor>>(30, s[i]);
    int iter = 0;
    for ( ; iter != info_->maxiter(); ++iter) {
      rall_[i]->zero();

      const double norm = t[i]->norm();
      t[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
        for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
          if (info_->sssr() && (jst != i || ist != i))
            continue;
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t[i]->at(ist);
          r = rall_[i]->at(jst);

          e0_ = e0all_[i] - info_->shift();
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
        }
      }
      // solve using subspace updates
      rall_[i] = solver->compute_residual(t[i], rall_[i]);
      t[i] = solver->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error = rall_[i]->rms();
      print_iteration(iter, energy_[i], error, mtimer.tick());
      conv = error < info_->thresh();

      // compute delta t2 and update amplitude
      if (!conv) {
        t[i]->zero();
        update_amplitude(t[i], rall_[i]);
      }
      if (conv) break;
    }
    if (i+1 != nstates_) cout << endl;
    converged &= conv;
  }
  print_iteration(!converged);
  return t;
}


void CASPT2::CASPT2::solve_deriv() {
  Timer timer;
  // First solve lambda equation if this is MS-CASPT2
  if (info_->do_ms() && nstates_ > 1) {
    // allocate lall_
    for (int i = 0; i != nstates_; ++i)
      lall_.push_back(t2all_[i]->clone());
    // lamda eqn :   T_M <omega' | H | M' > T_M' + <omega' | f - E0_M + Eshift | Omega> lamdba  - E_shift * (T_M)^2 * <proj|Psi_M>= 0
    // compute first term and shift term (if used)
    print_iteration();

    // rall_[0] stores the result of summation over M'
    const int target = info_->target();
    rall_[0]->zero();
    for (int istate = 0; istate != nstates_; ++istate) //K states
      rall_[0]->ax_plus_y((*heff_)(istate, target), sall_[istate]);

    for (int istate = 0; istate != nstates_; ++istate) { //K states
      sall_[istate]->zero();
      for (int jst = 0; jst != nstates_; ++jst)
        if (!info_->sssr() || istate == jst)
          sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, target), rall_[0]->at(jst));
      if (info_->shift() != 0.0) {
        if (!info_->shift_diag())
          throw runtime_error("Sorry, off-diagonal levelshift corrections are implemented for energies but not for gradients.");
        // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
        n = init_residual();
        for (int jst = 0; jst != nstates_; ++jst) { // bra
          for (int ist = 0; ist != nstates_; ++ist) { // ket
            if (info_->sssr() && (jst != istate || ist != istate))
              continue;
            set_rdm(jst, ist);
            t2 = t2all_[istate]->at(ist);
            shared_ptr<Queue> normq = make_normq(true, jst == ist);
            while (!normq->done())
              normq->next_compute();
            sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * pow((*heff_)(istate, target), 2.0), n);
          }
        }
      }
    }
    // solve linear equation and store lambda in lall
    lall_ = solve_linear(sall_, lall_);
    timer.tick_print("CASPT2 lambda equation");
  }

  if (lall_.empty()) {
    t2 = t2all_[0]->at(0);
    {
      den2 = h1_->clone();
      shared_ptr<Queue> dens2 = make_densityq();
      while (!dens2->done())
        dens2->next_compute();
      den2_ = den2->matrix();
    } {
      den1 = h1_->clone();
      shared_ptr<Queue> dens1 = make_density1q();
      while (!dens1->done())
        dens1->next_compute();
      den1_ = den1->matrix();
    } {
      Den1 = init_residual();
      shared_ptr<Queue> Dens1 = make_density2q();
      while (!Dens1->done())
        Dens1->next_compute();
      Den1_ = Den1;
    }
    timer.tick_print("Correlated density matrix evaluation");

    // first make ci_deriv_
    ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), 1);
    const size_t cisize = ci_deriv_->data(0)->size();

    const size_t cimax = 1000; // TODO interface to the input
    const size_t npass = (cisize-1)/cimax+1;
    const size_t chunk = (cisize-1)/npass+1;

    for (int ipass = 0; ipass != npass; ++ipass) {
      const size_t offset = ipass * chunk;
      const size_t size = min(chunk, cisize-offset);

      feed_rdm_deriv(offset, size); // this set ci_
      deci = make_shared<Tensor>(vector<IndexRange>{ci_});
      deci->allocate();
      shared_ptr<Queue> dec = make_deciq();
      while (!dec->done())
        dec->next_compute();
      copy_n(deci->vectorb()->data(), size, ci_deriv_->data(0)->data()+offset);

      stringstream ss; ss << "CI derivative evaluation " << setw(5) << ipass+1 << " / " << npass;
      timer.tick_print(ss.str());
    }
    cout << endl;
  } else {
    // in case when CASPT2 is not variational...
    MSCASPT2::MSCASPT2 ms(*this);
    ms.solve_deriv();
    den1_ = ms.rdm11();
    den2_ = ms.rdm12();
    Den1_ = ms.rdm21();
    ci_deriv_ = ms.ci_deriv();
    dcheck_ = ms.dcheck();
    timer.tick();
  }

  correlated_norm_.resize(nstates_);
  if (nstates_ == 1) {
    n = init_residual();
    shared_ptr<Queue> normq = make_normq();
    while (!normq->done())
      normq->next_compute();
    correlated_norm_[0] = dot_product_transpose(n, t2);
  } else {
    n = init_residual();
    for (int istate = 0; istate != nstates_; ++istate) {
      double tmp = 0.0;
      for (int jst = 0; jst != nstates_; ++jst) { // bra
        for (int ist = 0; ist != nstates_; ++ist) { // ket
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          shared_ptr<Queue> normq = make_normq(true, jst == ist);
          while (!normq->done())
            normq->next_compute();
          tmp += dot_product_transpose(n, lall_[istate]->at(jst));
        }
      }
      correlated_norm_[istate] = tmp;
    }
  }
  timer.tick_print("T1 norm evaluation");

  // some additional terms to be added
  const int ncore = info_->ncore();
  const int nclosed = info_->nclosed()-info_->ncore();
  const int nact = info_->nact();
  {
    // d_1^(2) -= <1|1><0|E_mn|0>     [Celani-Werner Eq. (A6)]
    auto dtmp = den2_->copy();
    for (int ist = 0; ist != nstates_; ++ist) {
      auto rdmtmp = rdm1all_->at(ist, ist)->matrix();
      for (int i = nclosed; i != nclosed+nact; ++i)
        for (int j = nclosed; j != nclosed+nact; ++j)
          dtmp->element(j, i) -=  correlated_norm_[ist] * (*rdmtmp)(j-nclosed, i-nclosed);
    }
    dtmp->symmetrize();
    den2_ = dtmp;
  }

  shared_ptr<const Reference> ref = info_->ref();
  const MatView acoeff = coeff_->slice(nclosed+ncore, nclosed+ncore+nact);

  // code to calculate h+g(d). When add is false, h is not added
  auto focksub = [&](shared_ptr<const Matrix> moden, const MatView coeff, const bool add) {
    shared_ptr<const Matrix> jop = ref->geom()->df()->compute_Jop(make_shared<Matrix>(coeff * *moden ^ coeff));
    auto out = make_shared<Matrix>(acoeff % (add ? (*ref->hcore() + *jop) : *jop) * acoeff);
    shared_ptr<const DFFullDist> full = ref->geom()->df()->compute_half_transform(acoeff)->compute_second_transform(coeff)->apply_J()->swap();
    shared_ptr<DFFullDist> full2 = full->copy();
    full2->rotate_occ1(moden);
    *out += *full->form_2index(full2, -0.5);
    return out;
  };
  shared_ptr<const Matrix> fock = focksub(ref->rdm1_mat(), coeff_->slice(0, ref->nocc()), true); // f
  {
    // correct cideriv for fock derivative [Celani-Werner Eq. (C1), some terms in first and second lines]
    // y_I += (g[d^(2)]_ij - Nf_ij) <I|E_ij|0>
    shared_ptr<const Matrix> gd2 = focksub(den2_, coeff_->slice(ncore, coeff_->mdim()), false); // g(d2)

    for (int ist = 0; ist != nstates_; ++ist) {
      const Matrix op(*gd2 * (1.0/nstates_) - *fock * correlated_norm_[ist]);
      shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          ci_deriv_->data(ist)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
    }

    // y_I += <I|H|0> (for mixed states); taking advantage of the fact that unrotated CI vectors are eigenvectors
    const int tst = info_->target();
    const Matrix ur(xmsmat_ ? *xmsmat_ * *heff_ : *heff_);
    for (int ist = 0; ist != nstates_; ++ist)
      for (int jst = 0; jst != nstates_; ++jst)
        ci_deriv_->data(jst)->ax_plus_y(2.0*ur(ist,tst)*(*heff_)(jst,tst)*ref->energy(ist), info_orig_->ciwfn()->civectors()->data(ist));
  }

  // finally if this is XMS-CASPT2 gradient computation, we compute dcheck and contribution to y
  if (xmsmat_) {
    Matrix wmn(nstates_, nstates_);
    shared_ptr<Tensor> dc = rdm1_->clone();
    for (int i = 0; i != nstates_; ++i)
      for (int j = 0; j != i; ++j) {
        const double cy = info_->ciwfn()->civectors()->data(j)->dot_product(ci_deriv_->data(i))
                        - info_->ciwfn()->civectors()->data(i)->dot_product(ci_deriv_->data(j));
        wmn(j,i) = fabs(e0all_[j]-e0all_[i]) > 1.0e-12 ? -0.5 * cy / (e0all_[j]-e0all_[i]) : 0.0;
        wmn(i,j) = wmn(j,i);
        dc->ax_plus_y(wmn(j,i), rdm1all_->at(j, i));
        dc->ax_plus_y(wmn(i,j), rdm1all_->at(i, j));
      }
    dcheck_ = dc->matrix();

    // fill this into CI derivative. (Y contribution is done inside Z-CASSCF together with frozen core)
    shared_ptr<const Matrix> gdc = focksub(dcheck_, acoeff, false);
    for (int ist = 0; ist != nstates_; ++ist) {
      shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
      for (int jst = 0; jst != nstates_; ++jst) {
        Matrix op(*fock * wmn(jst, ist));
        if (ist == jst)
          op += *gdc * (1.0/nstates_);
        for (int i = 0; i != nact; ++i)
          for (int j = 0; j != nact; ++j)
            ci_deriv_->data(jst)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
      }
    }

    // also rotate cideriv back to the MS states
    btas::contract(1.0, *ci_deriv_->copy(), {0,1,2}, (*xmsmat_), {3,2}, 0.0, *ci_deriv_, {0,1,3});
  }

  // restore original energy
  energy_ = pt2energy_;
  timer.tick_print("Postprocessing SMITH");
}

#endif
