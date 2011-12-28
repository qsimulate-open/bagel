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
#include <src/scf/scf.h>
#include <src/util/f77.h>

using namespace std;

SuperCI::SuperCI(const multimap<string, string> idat, const shared_ptr<Geometry> geom)
 : CASSCF(idat, geom) {
  common_init();
}

SuperCI::SuperCI(const multimap<string, string> idat, const shared_ptr<Geometry> geom, const shared_ptr<SCF> scf)
 : CASSCF(idat, geom, scf) {
  common_init();
}

void SuperCI::common_init() {
  cout << "    * using the Super CI algorithm as noted in Roos (1980) IJQC" << endl << endl;
}

SuperCI::~SuperCI() {

}


void SuperCI::compute() {
  const string indent = "  ";

  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ")===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  // initializing Hcore matrix (redundant copy, but I can live with it).
  shared_ptr<Fock> hcore_;
  {
    shared_ptr<Hcore> hc(new Hcore(geom_));
    shared_ptr<Fock> fc(new Fock(geom_, hc)); hcore_ = fc;
  }
  shared_ptr<RotFile> qxr(new RotFile(nclosed_, nact_, nvirt_)); 

  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();

    // first perform CASCI to obtain RDMs
    mute_stdcout();
    fci_->compute();
    fci_->compute_rdm12();
    resume_stdcout();
    // get energy
    std::vector<double> energy = fci_->energy();
    // slot in density matrix
    shared_ptr<Matrix1e> denall = ao_rdm1(fci_->rdm1_av());
    shared_ptr<Matrix1e> deninact = ao_rdm1(fci_->rdm1_av(), true); // true means inactive_only

    // here make a natural orbitals and update the coefficients
    // this effectively updates 1,2RDM and integrals 
    pair<vector<double>, vector<double> > natorb = fci_->natorb_convert();
    // new coefficients
    ref_->set_coeff(update_coeff(ref_->coeff(), natorb.first));
    // occupation number of the natural orbitals
    occup_ = natorb.second;

    // get quantity Q_xr = (xs|tu)P_rs,tu (x=general)
    // note: this should be after natorb transformation.
// strictly speaking, this line is not needed TODO
    qxr->zero();
    compute_qxr(fci_->jop()->mo2e_1ext_ptr(), fci_->rdm2_av(), qxr);

    shared_ptr<RotFile> cc_(new RotFile(nclosed_, nact_, nvirt_));
    {
      // Davidson utility. We diagonalize a super CI matrix every macro iteration
      DavidsonDiag<RotFile> davidson(1, max_iter_);
      shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_));

      // first, <proj|H|0> is computed 
      sigma_->zero();
      cc_->zero();
      cc_->ele_ref() = 1.0;

      // computes f and f_act
      // TODO call to fock builder should be done only once.
      shared_ptr<Matrix1e> f, finact;
      shared_ptr<Coeff> coeff = ref_->coeff();
      {
        shared_ptr<Fock> f_ao(new Fock(geom_, hcore_, denall, ref_->schwarz()));
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *f_ao * *coeff)); f = ft;
      } {
        shared_ptr<Fock> finact_ao(new Fock(geom_, hcore_, deninact, ref_->schwarz()));
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *finact_ao * *coeff)); finact = ft;
      }

      // <a/i|H|0> = 2f_ai
      grad_vc(f, sigma_);
      // <a/r|H|0> = f^inact_as d_sr + 2(as|tu)P_rs,tu 
      grad_va(finact, fci_->rdm1_av(), qxr, sigma_);
      // <r/i|H|0> = 2f_ri - f^inact_is d_sr - 2(is|tu)P_rs,tu
      grad_ca(f, finact, fci_->rdm1_av(), qxr, sigma_);

      // then microiteration for diagonalization 
      for (int miter = 0; miter != max_micro_iter_; ++miter) {

        if (miter != 0) {
#if 0
          // equation 21d
          sigma_ai_ai_;

          // equation 21e
          sigma_ai_at_;
          sigma_at_ai_;

          // equation 21f
          sigma_at_at_;

          // equation 21b
          sigma_ai_ti_;
          sigma_ti_ai_;

          // equation 21a
          sigma_ti_ti_;
#endif
        }

        shared_ptr<RotFile> ccn(new RotFile(cc_));
        vector<shared_ptr<RotFile> > ccp(1,ccn);
        shared_ptr<RotFile> sigman(new RotFile(sigma_));
        vector<shared_ptr<RotFile> > sigmap(1,sigman);
        const vector<double> energies = davidson.compute(ccp, sigmap);
        const double energy = energies.front();
        cout << "   state averaged energy " << energy << endl;
      }
    }

    // compute errors
    std::vector<double> error(nstate_, 0.0);
    for (int i = 0; i != nstate_; ++i) {
//TODO errors.push_back(errvec[i]->variance());
      conv[i] = static_cast<int>(error[i] < thresh_);
    }

    int end = ::clock();
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      cout << indent << setw(5) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                                        << setw(17) << fixed << setprecision(8) << energy[i] << "   " 
                                        << setw(10) << scientific << setprecision(2) << error[i] << fixed << setw(10) << setprecision(2)
                                        << (end - start)/static_cast<double>(CLOCKS_PER_SEC) << endl; 
    }

#if 0
    if (*min_element(conv.begin(), conv.end())) break;
#endif
    if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in the CASSCF macro interation." << endl << endl;
      break;
    }
  }

}


shared_ptr<Coeff> SuperCI::update_coeff(const shared_ptr<Coeff> cold, vector<double> mat) {
  shared_ptr<Matrix1e> cnew(new Matrix1e(*dynamic_cast<Matrix1e*>(cold.get())));
  int nbas = geom_->nbasis();
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, &(mat[0]), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  shared_ptr<Coeff> cnew1(new Coeff(*cnew)); 
  return cnew1;
}


// <a/i|H|0> = 2f_ai
void SuperCI::grad_vc(const shared_ptr<Matrix1e> fock, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  fill(target, target+nclosed_*nvirt_, 0.0);
  for (int i = 0; i != nclosed_; ++i, target += nvirt_)
    daxpy_(nvirt_, 2.0, fock->element_ptr(nclosed_+nact_,i), 1, target, 1); 
}


// <a/r|H|0> = (f-f^act)_as d_sr + 2(as|tu)P_rs,tu 
void SuperCI::grad_va(const shared_ptr<Matrix1e> finact, shared_ptr<RDM<1> > den,
                      shared_ptr<RotFile> qxr, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nact_) return;
  copy(qxr->ptr_va(), qxr->ptr_va()+nvirt_*nact_, sigma->ptr_va());
  dgemm_("N", "N", nvirt_, nact_, nact_, 1.0, finact->data()+nclosed_*nbasis_+nclosed_+nact_, nbasis_, den->first(), nact_,
                   1.0, sigma->ptr_va(), nvirt_);
}


// <r/i|H|0> = 2f_ri - f^inact_is d_sr - 2(is|tu)P_rs,tu
void SuperCI::grad_ca(const shared_ptr<Matrix1e> f, const shared_ptr<Matrix1e> finact, shared_ptr<RDM<1> > den,
                      shared_ptr<RotFile> qxr, shared_ptr<RotFile> sigma) {
  if (!nclosed_ || !nact_) return;
  const int n = nclosed_*nact_;
  double* target = sigma->ptr_ca();
  fill(target, target+n, 0.0);
  for (int i = 0; i != nact_; ++i, target += nclosed_)
    daxpy_(nclosed_, 2.0, f->element_ptr(0,nclosed_+i), 1, target, 1);
  dgemm_("N", "N", nclosed_, nact_, nact_, -1.0, finact->data()+nclosed_*nbasis_, nbasis_, den->first(), nact_,
                   1.0, sigma->ptr_ca(), nclosed_); 
  daxpy_(n, -1.0, qxr->ptr_ca(), 1, sigma->ptr_ca(), 1); 
}


void SuperCI::compute_qxr(double* int1ext, shared_ptr<RDM<2> > rdm2, shared_ptr<RotFile> qxr) {
  // int1ext = (st|ux), rdm2 = D_st,ur
  const int nbas = geom_->nbasis(); // caution :: this is AO and therefore not nbasis_
  const int common = nact_*nact_*nact_; 
  double* buf = new double[nbas*nact_];
  dgemm_("T", "N", nbas, nact_, common, 1.0, int1ext, common, rdm2->first(), common, 0.0, buf, nbas); 
  // slot in to apropriate places
  if (nclosed_)
  dgemm_("T", "N", nclosed_, nact_, nbas, 1.0, ref_->coeff()->data(), nbas, buf, nbas, 0.0, qxr->ptr_ca(), nclosed_); 
  if (nvirt_) 
  dgemm_("T", "N", nvirt_,   nact_, nbas, 1.0, ref_->coeff()->data()+nbas*(nclosed_+nact_),
                                                                      nbas, buf, nbas, 0.0, qxr->ptr_va(), nvirt_); 
  delete[] buf;
}
