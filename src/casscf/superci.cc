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
#include <src/casscf/f77.h>

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

void SuperCI::update_orbitals_() {

}

void SuperCI::update_civectors_() {

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
  shared_ptr<Coeff> coeff = ref_->coeff();
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
    shared_ptr<Matrix1e> denact = ao_rdm1(fci_->rdm1_av(), true); // true means active_only

    // get quantity Q_xr = (xs|tu)P_rs,tu (x=general)
    qxr->zero();
    compute_qxr(fci_->jop()->mo2e_1ext_ptr(), fci_->rdm2_av(), qxr);
    qxr->print();

    shared_ptr<RotFile> cc_(new RotFile(nclosed_, nact_, nvirt_));
    {
      // Davidson utility. We diagonalize a super CI matrix every macro iteration
      DavidsonDiag<RotFile> davidson(nstate_, max_iter_);
      shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_));

      // first, <proj|H|0> is computed 
      cc_->ele_ref() = 1.0;
      sigma_->ele_ref() = 0.0;

      // computes f and f_act
      // TODO call to fock builder should be only once.
      shared_ptr<Matrix1e> f, fact;
      {
        shared_ptr<Fock> f_ao(new Fock(geom_, hcore_, denall, ref_->shwarz()));
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *f_ao * *coeff)); f = ft;
      } {
        shared_ptr<Fock> fact_ao(new Fock(geom_, hcore_, denact, ref_->shwarz()));
        shared_ptr<Matrix1e> ft(new Matrix1e(*coeff % *fact_ao * *coeff)); fact = ft;
      }

      // <a/i|H|0> = 2f_ai
      grad_vc(f, sigma_);
      // <a/r|H|0> = (f-f^act)_as d_sr + 2(as|tu)P_rs,tu 
      double* data_va = sigma_->ptr_va();
//    grad_va(f, fact, sigma_);
      // <r/i|H|0> = 2f_ri - (f-f^act)_is d_sr - 2(is|tu)P_rs,tu
      double* data_ca = sigma_->ptr_ca();

      // then microiteration for diagonalization 
      for (int miter = 0; miter != max_micro_iter_; ++miter) {
        //
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


// <a/i|H|0> = 2f_ai
void SuperCI::grad_vc(const shared_ptr<Matrix1e> fock, shared_ptr<RotFile> sigma) {
  double* target = sigma->ptr_vc(); // v runs first
  const double two = 2.0;
  const int unit = 1;
  for (int i = 0; i != nclosed_; ++i, target += nvirt_)
    daxpy_(&nvirt_, &two, fock->element_ptr(nclosed_+nact_,i), &unit, target, &unit); 
}


void SuperCI::compute_qxr(double* int1ext, shared_ptr<RDM<2> > rdm2, shared_ptr<RotFile> qxr) {
  // int1ext = (st|ux), rdm2 = D_st,ur , !! factor of 0.5 included here
  const int nbas = geom_->nbasis(); // caution :: this is AO and therefore not nbasis_
  const int common = nact_*nact_*nact_; 
  const double one = 1.0;
  const double half = 0.5;
  const double zero = 0.0; 
  double* buf = new double[nbas*nact_];
  dgemm_("T", "N", &nbas, &nact_, &common, &half, int1ext, &common, rdm2->first(), &common, &zero, buf, &nbas); 
  // slot in to apropriate places
  dgemm_("T", "N", &nclosed_, &nact_, &nbas, &one, ref_->coeff()->data(), &nbas, buf, &nbas, &zero, qxr->ptr_ca(), &nclosed_); 
  dgemm_("T", "N", &nvirt_,   &nact_, &nbas, &one, ref_->coeff()->data()+nbas*(nclosed_+nact_),
                                                                          &nbas, buf, &nbas, &zero, qxr->ptr_va(), &nvirt_); 
  delete[] buf;
}
