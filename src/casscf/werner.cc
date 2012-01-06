//
// Author : Toru Shiozaki
// Date   : Jan 2012
//


#include <src/casscf/werner.h>
#include <src/casscf/jvec.h>

using namespace std;

#define DF 1

void WernerKnowles::compute() {
  const string indent = "  ";

  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

  // initializing Hcore matrix (redundant copy, but I can live with it).
  shared_ptr<Fock<DF> > hcore_;
  {
    shared_ptr<Hcore> hc(new Hcore(geom_));
    shared_ptr<Fock<DF> > fc(new Fock<DF>(geom_, hc)); hcore_ = fc;
  }

  // macro iteration
  for (int iter = 0; iter != 1; ++iter) {
//for (int iter = 0; iter != max_iter_; ++iter) {

    // performs FCI, which also computes one-index transformed integrals 
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    vector<double> energy = fci_->energy();

    shared_ptr<Jvec> jvec(new Jvec(fci_, ref_->coeff(), nclosed_, nact_, nvirt_));

    for (int miter = 0; miter != max_micro_iter_; ++miter) {

    }

  }

}

