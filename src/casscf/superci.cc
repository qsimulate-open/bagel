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

    shared_ptr<RotFile> cc_(new RotFile(nclosed_, nact_, nvirt_));
    {
      // Davidson utility. We diagonalize a super CI matrix every macro iteration
      DavidsonDiag<RotFile> davidson(nstate_, max_iter_);
      shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_));

      // first, <proj|H|0> is computed 
      cc_->ele_ref() = 1.0;
      sigma_->ele_ref() = 0.0;
      // will compute f and f_act
      shared_ptr<Fock> f(new Fock(geom_, hcore_, denall, ref_->shwarz()));
      shared_ptr<Fock> fact(new Fock(geom_, hcore_, denact, ref_->shwarz()));
      // <a/i|H|0> = 2f_ai
      double* data_vc = sigma_->ptr_vc();
      // <a/r|H|0> = (f-f^act)_as d_sr + 2(as|tu)P_rs,tu 
      double* data_va = sigma_->ptr_va();
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

    if (*min_element(conv.begin(), conv.end())) break;
    if (iter == max_iter_-1) {
      cout << indent << endl << indent << "  * Max iteration reached in the CASSCF macro interation." << endl << endl;
      break;
    }
  }

}
