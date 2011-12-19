//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <iostream>
#include <stdexcept>
#include <src/fci/mofile.h>
#include <src/fci/fci.h>
#include <src/fci/civec.h>

using namespace std;

static const int unit = 1;

void FCI::compute() {

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment."); 

  print_header();

  // first obtain reference function
  ref_->compute();
  shared_ptr<Coeff> cmo = ref_->coeff(); 

  // iiii file to be created (MO transformation).
  // now Jop->mo1e() and Jop->mo2e() contains one and two body part of Hamiltonian
  shared_ptr<MOFile> Jop(new MOFile(geom_, cmo));

  // right now full basis is used. 
  Jop->create_Jiiii(ncore_, ncore_+norb_);

  // hcore weighted coupling lists

  // Creating an initial CI vector
  shared_ptr<Civec> cc(new Civec(stringb_.size(), stringa_.size())); // B runs first
  shared_ptr<Civec> sigma(new Civec(stringb_.size(), stringa_.size()));

  // TODO this is only if HF is valid and therefore wrong in general cases.
  cc->element(0,0) = 1.0;

  // main iteration starts here

  // one-electron alpha: sigma(Psi'a, Psib) += sign h'(ij) C(Psia, Psib) 
#if 0
  {
    const int lb = stringb_.size(); 
    for (auto iter = phia_.begin();  iter != phia_.end(); ++iter) {
      const double hc = Jop->mo1e(get<2>(*iter)) * get<1>(*iter);
      daxpy_(&lb, &hc, &(cc->element(0, get<3>(*iter))), &unit, &(sigma->element(0, get<0>(*iter))), &unit); 
    }
  }
  // step (c)

  // transpose the coefficient matrix
  shared_ptr<Civec> ct = cc->transpose();


  // transpose the sigma matrix

  // one-electron beta: sigma(Psia, Psib') += sign h'(ij) C(Psia, Psib) 
  {
    const int lb = stringa_.size(); 
    for (auto iter = phib_.begin();  iter != phib_.end(); ++iter) {
      const double hc = Jop->mo1e(get<2>(*iter)) * get<1>(*iter);
      daxpy_(&lb, &hc, &(cc->element(0, get<3>(*iter))), &unit, &(sigma->element(0, get<0>(*iter))), &unit); 
    }
  }
#endif


  // main iteration ends here
}

void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}
