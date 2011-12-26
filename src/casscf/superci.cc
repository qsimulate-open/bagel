//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/superci.h>
#include <iostream>
#include <src/fci/fci.h>
#include <src/casscf/rotfile.h>

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

  cout << "    * Initial CASCI RDMs are generated" << endl << endl;
  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ")===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  shared_ptr<RotFile> cc_(new RotFile(nclosed_, nact_, nvirt_));
  shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_));

  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();

    // first perform CASCI to obtain RDMs
    mute_stdcout();
    fci_->compute();
    fci_->compute_rdm12();
    resume_stdcout();
    // get energy
    std::vector<double> energy = fci_->energy();

    int nmicro = 10;

    // first, 

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
