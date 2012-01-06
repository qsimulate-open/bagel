//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/superci.h>
#include <iostream>
#include <src/fci/fci.h>
#include <src/casscf/rotfile.h>
#include <src/util/davidson.h>
#include <src/scf/hcore.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/util/hpw_diis.h>

using namespace std;

#define DF 1

static const double cps = static_cast<double>(CLOCKS_PER_SEC);

void SuperCI::compute() {
  const string indent = "  ";

  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

  // initializing Hcore matrix (redundant copy, but I can live with it).
  shared_ptr<Fock<DF> > hcore_;
  {
    shared_ptr<Hcore> hc(new Hcore(geom_));
    shared_ptr<Fock<DF> > fc(new Fock<DF>(geom_, hc)); hcore_ = fc;
  }

  // DIIS: will be turned on at iter = diis_start_ (>1), 
  //       update log(U) where Cnow = Corig U. This is basically the same as the Hampel-Peterson-Werner
  //       paper on Brueckner CC
  shared_ptr<HPW_DIIS<Matrix1e> > diis;

  // ============================
  // macro iteration from here
  // ============================
  double gradient = 1.0e100;
  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();

    if (iter >= diis_start_ && gradient < 1.0e-4 && !diis) {
      shared_ptr<Matrix1e> tmp(new Matrix1e(*ref_->coeff()));
      shared_ptr<HPW_DIIS<Matrix1e> > d(new HPW_DIIS<Matrix1e>(5, tmp));
      diis = d;
    }

    // first perform CASCI to obtain RDMs
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    vector<double> energy = fci_->energy();

    int start0 = ::clock();

    // here make a natural orbitals and update the coefficients
    // this effectively updates 1,2RDM and integrals
    const pair<vector<double>, vector<double> > natorb = fci_->natorb_convert();
    // new coefficients
    shared_ptr<Coeff> new_coeff = update_coeff(ref_->coeff(), natorb.first);
    ref_->set_coeff(new_coeff);
    // occupation number of the natural orbitals
    occup_ = natorb.second;
    if (std::abs(occup_.front()-2.0) < 1.0e-16 || std::abs(occup_.back()) < 1.0e-16)
      throw runtime_error("CASSCF does not work so far if occupied orbitals are strictly doubly occupied or empty.");

    // get quantity Q_xr = 2(xs|tu)P_rs,tu (x=general)
    // note: this should be after natorb transformation.
    shared_ptr<Qvec> qxr(new Qvec(geom_->nbasis(), nact_, geom_->df(), ref_->coeff(), nclosed_, fci_));

    cout << "     * Natural orbital transformation + Q vector  " << setprecision(2) << (::clock() - start0)/cps << " sec" << endl; start0 = ::clock();

    shared_ptr<RotFile> cc_(new RotFile(nclosed_, nact_, nvirt_));
    
    // Davidson utility. We diagonalize a super CI matrix every macro iteration
    DavidsonDiag<RotFile> davidson(1, max_micro_iter_);
    shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_));

    // computes f and f_act
    // TODO call to fock builder should be done only once.
    shared_ptr<Matrix1e> f, finact;
    shared_ptr<QFile>    fact, factp, gaa;
    shared_ptr<Coeff> coeff = ref_->coeff();

    {
      // Fock operators
      if (nclosed_) {
        shared_ptr<Matrix1e> deninact = ao_rdm1(fci_->rdm1_av(), true); // true means inactive_only
        shared_ptr<Matrix1e> f_inactao(new Matrix1e(geom_));
        dcopy_(f_inactao->size(), fci_->jop()->core_fock_ptr(), 1, f_inactao->data(), 1);
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *f_inactao * *coeff));
        finact = ft;

        shared_ptr<Matrix1e> denall = ao_rdm1(fci_->rdm1_av());
        shared_ptr<Matrix1e> denact(new Matrix1e(*denall-*deninact));
        shared_ptr<Fock<DF> > fact_ao(new Fock<DF>(geom_, hcore_, denact, ref_->schwarz()));
        shared_ptr<Matrix1e> fact(new Matrix1e(*coeff % (*fact_ao-*hcore_) * *coeff));
        shared_ptr<Matrix1e> fs(new Matrix1e(*finact+*fact));
        f = fs;
      } else {
        shared_ptr<Matrix1e> denall = ao_rdm1(fci_->rdm1_av());
        shared_ptr<Fock<DF> > f_ao(new Fock<DF>(geom_, hcore_, denall, ref_->schwarz()));
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *f_ao * *coeff)); f = ft;

        finact = hcore_;
      }
    }
    cout << "     * Fock + inactive Fock                       " << setprecision(2) << (::clock() - start0)/cps << " sec" << endl; start0 = ::clock();
    {
      // active-x Fock operator Dts finact_sx + Qtx
      shared_ptr<QFile> ft(new QFile(*qxr)); fact = ft; // nbasis_ runs first
      if (nclosed_) {
        for (int i = 0; i != nact_; ++i)
          daxpy_(nbasis_, occup_[i], finact->element_ptr(0,nclosed_+i), 1, fact->data()+i*nbasis_, 1);
      }
    }
    {
      // active Fock' operator (Fts+Fst) / (ns+nt)
      shared_ptr<QFile> ft(new QFile(nact_, nact_)); factp = ft;
      for (int i = 0; i != nact_; ++i)
        for (int j = 0; j != nact_; ++j)
          factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) / (occup_[i]+occup_[j]);
    }

    // G matrix (active-active) 2Drs,tu Factp_tu - delta_rs nr sum_v Factp_vv
    shared_ptr<QFile> ft4(new QFile(nact_, nact_)); gaa = ft4;
    dgemv_("N", nact_*nact_, nact_*nact_, 1.0, fci_->rdm2_av()->data(), nact_*nact_, factp->data(), 1, 0.0, gaa->data(), 1);
    double p = 0.0;
    for (int i = 0; i != nact_; ++i) p += occup_[i] * factp->element(i,i);
    for (int i = 0; i != nact_; ++i) gaa->element(i,i) -= occup_[i] * p;

    // first, <proj|H|0> is computed
    sigma_->zero();
    cc_->zero();
    cc_->ele_ref() = 1.0;

    // <a/i|H|0> = 2f_ai
    grad_vc(f, sigma_);
    // <a/r|H|0> = h_as d_sr + 2(as|tu)P_rs,tu = fact_rs
    grad_va(fact, sigma_);
    // <r/i|H|0> = 2f_ri - f^inact_is d_sr - 2(is|tu)P_rs,tu = 2f_ri - fact_ri
    grad_ca(f, fact, sigma_);
    sigma_->ele_ref() = 0.0;

    // setting error of macro iteration
    gradient = sigma_->ddot(*sigma_) / sigma_->size();

    // denominator
    shared_ptr<RotFile> denom_ = const_denom(gaa, factp, f);

    shared_ptr<RotFile> init_sigma(new RotFile(*sigma_));

    // ---------------------------------------
    // then microiteration for diagonalization
    // ---------------------------------------
    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      const int mstart = ::clock();

      if (miter != 0) {
        sigma_->zero();

        // equation 21d
        sigma_ai_ai_(cc_, sigma_, f);
        // equation 21e
        sigma_at_ai_(cc_, sigma_, fact);
        // equation 21f // note a typo!
        sigma_at_at_(cc_, sigma_, gaa, f);
        // equation 21b
        sigma_ai_ti_(cc_, sigma_, fact);
        // equation 21a
        sigma_ti_ti_(cc_, sigma_, gaa, f, factp);

        // projection to reference
        cc_->ele_ref()=0.0;
        sigma_->ele_ref() = init_sigma->ddot(*cc_);
      }

      // enters davidson iteration
      shared_ptr<RotFile> ccp(new RotFile(*cc_));
      shared_ptr<RotFile> sigmap(new RotFile(*sigma_));
      const double mic_energy = davidson.compute(ccp, sigmap);

      // residual vector and error
      shared_ptr<RotFile> residual = davidson.residual().front();
      const double error = residual->ddot(*residual) / residual->size();

      const int mend = ::clock();
      if (miter == 0) cout << endl << "     == micro iteration == " << endl;
      cout << setw(10) << miter << "   " << setw(20) << setprecision(12) << mic_energy << " "
           << setw(10) << scientific << setprecision(2) << error << fixed << " " << (mend - mstart)/cps << endl;

      if (error < thresh_micro_) { cout << endl; break; }
      if (miter+1 == max_micro_iter_) throw runtime_error("max_micro_iter_ is reached in CASSCF");


      // update cc_
      for (double *i = residual->begin(), *j = denom_->begin(); i != residual->end(); ++i, ++j) { *i /= *j; }
      const double a = davidson.orthog(residual);
      cc_ = residual;
    }
    // ---------------------------------------
    // micro iteration to here
    // ---------------------------------------


    // rotation parameters
    cc_ = davidson.civec().front();
    dscal_(cc_->size()-1, 1.0/cc_->ele_ref(), cc_->data(), 1);
    // unitary matrix
    shared_ptr<Matrix1e> rot = cc_->unpack(ref_->coeff()->geom())->exp();
    // forcing rot to be unitary (usually not needed, though)
    rot->purify_unitary();

    if (!diis) {
      *ref_->coeff() *= *rot;
    } else {
      // including natorb.first to rot so that they can be processed at once
      shared_ptr<Matrix1e> tmp(new Matrix1e(*rot));
      dgemm_("N", "N", nact_, nbasis_, nact_, 1.0, &(natorb.first[0]), nact_, rot->element_ptr(nclosed_, 0), nbasis_, 0.0,
                                                                       tmp->element_ptr(nclosed_, 0), nbasis_);
      tmp = tailor_rotation(tmp);
      
      shared_ptr<Matrix1e> mcc = diis->extrapolate(tmp);
      shared_ptr<Coeff> newcc(new Coeff(*mcc));
      ref_->set_coeff(newcc);
    }

    // print out...
    int end = ::clock();
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      resume_stdcout();
      cout << indent << setw(5) << iter << setw(3) << i << setw(2) << (diis ? "*" : " ") 
                                        << setw(17) << fixed << setprecision(8) << energy[i] << "   "
                                        << setw(10) << scientific << setprecision(2) << (i==0 ? gradient : 0.0) << fixed << setw(10) << setprecision(2)
                                        << (end - start)/cps << endl;
      mute_stdcout();
    }

    if (gradient < thresh_) break;
    if (iter == max_iter_-1) {
      resume_stdcout();
      cout << indent << endl << indent << "  * Max iteration reached in the CASSCF macro interation." << endl << endl;
      mute_stdcout();
      break;
    }
  }
  // ============================
  // macro iteration to here
  // ============================
  resume_stdcout();

}


// rotate (within allowed rotations) the transformation matrix so that it is diagonal in each subblock 
shared_ptr<Matrix1e> SuperCI::tailor_rotation(const shared_ptr<Matrix1e> seed) {

  shared_ptr<Matrix1e> out = seed->clone();
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nclosed_; ++j)
      out->element(j,i) = seed->element(j,i);
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      out->element(j+nclosed_,i+nclosed_) = seed->element(j+nclosed_,i+nclosed_);
  for (int i = 0; i != nvirt_; ++i)
    for (int j = 0; j != nvirt_; ++j)
      out->element(j+nocc_,i+nocc_) = seed->element(j+nocc_,i+nocc_);
  out->inverse(); 
  out->purify_unitary();
  *out = *seed * *out;

  return out;
}


shared_ptr<Coeff> SuperCI::update_coeff(const shared_ptr<Coeff> cold, vector<double> mat) const {
  shared_ptr<Matrix1e> cnew(new Matrix1e(*dynamic_cast<Matrix1e*>(cold.get())));
  int nbas = geom_->nbasis();
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, &(mat[0]), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  shared_ptr<Coeff> cnew1(new Coeff(*cnew));
  return cnew1;
}




void SuperCI::compute_qxr(double* int1ext, shared_ptr<RDM<2> > rdm2, shared_ptr<QFile> qxr) {
  // int1ext = (xu|st) = (ts|ux), rdm2 = D_ru,st = D_ur,ts = D_ts,ur
  const int nbas = geom_->nbasis(); // caution :: this is AO and therefore not nbasis_
  const int common = nact_*nact_*nact_;
  QFile buf(nbas,nact_);
  dgemm_("T", "N", nbas, nact_, common, 1.0, int1ext, common, rdm2->first(), common, 0.0, buf.data(), nbas);
  dgemm_("T", "N", nbasis_, nact_, nbas, 1.0, ref_->coeff()->data(), nbas, buf.data(), nbas, 0.0, qxr->data(), nbasis_);
}


// compute denominators
shared_ptr<RotFile> SuperCI::const_denom(const shared_ptr<QFile> gaa, const shared_ptr<QFile> factp, const shared_ptr<Matrix1e> f) const {
  shared_ptr<RotFile> denom(new RotFile(nclosed_, nact_, nvirt_));
  fill(denom->data(), denom->data()+denom->size(), 1.0e100);

  double* target = denom->ptr_va();
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_; ++j, ++target) {
      *target = gaa->element(i,i) / occup_[i] + f->element(j+nocc_, j+nocc_);
    }
  }

  target = denom->ptr_vc();
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_; ++j, ++target) {
      *target = f->element(j+nocc_, j+nocc_) - f->element(i, i);
    }
  }

  target = denom->ptr_ca();
  for (int i = 0; i != nact_; ++i) {
    const double fac = -((2.0 - 2.0*occup_[i]) * factp->element(i, i) - gaa->element(i, i)) / (2.0 - occup_[i]);
    for (int j = 0; j != nclosed_; ++j, ++target) {
      *target = fac - f->element(j, j);
    }
  }
  return denom;
}


