//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <src/fci/fci.h>
#include <src/util/davidson.h>

// TODO hardwired
#define MAX_ITER_FCI 100
#define THREASH 1.0e-8 // variant

using namespace std;

static const int unit = 1;
static const double one = 1.0;
static const double half = 0.5;
static const double zero = 0.0; 
static const string indent = "  ";
static const string space3 = "   "; 

void FCI::compute() {
  const int num_state = 1; // TODO should be read from the input

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment."); 

  // first obtain reference function
  ref_->compute();
  shared_ptr<Coeff> cmo = ref_->coeff(); 

  // iiii file to be created (MO transformation).
  // now Jop->mo1e() and Jop->mo2e() contains one and two body part of Hamiltonian
  shared_ptr<MOFile> Jop(new MOFile(geom_, cmo));

  // right now full basis is used. 
  Jop->create_Jiiii(ncore_, ncore_+norb_);

  // create denominator here. Stored in shared<Civec> denom_
  const_denom(Jop);

  // some constants
  const int la = stringa_.size();
  const int lb = stringb_.size();
  const int ij = norb_ * norb_; 

  // we need two vectors for intermediate quantities
  shared_ptr<Dvec> d(new Dvec(lb, la, ij));
  shared_ptr<Dvec> e(new Dvec(lb, la, ij));

  // Creating an initial CI vector
  shared_ptr<Dvec> cc(new Dvec(lb, la, num_state)); // B runs first
  shared_ptr<Dvec> sigma(new Dvec(lb, la, num_state));

  // TODO This is only if RHF is valid and therefore wrong in general cases.
  //      The right ways is to compute diagonal dinominators and select determinants
  //      that have small values, and then spin adapt.
  // TODO multiple state runs should be considered. At this moment, it is only partially..
  cc->data(0)->element(0,0) = 1.0;

  // nuclear energy retrieved from geometry
  const double nuclear = geom_->nuclear_repulsion();

  // Davidson utility
  DavidsonDiag<Civec> davidson(num_state, MAX_ITER_FCI);

  // main iteration starts here
  cout << "  === FCI iteration ===" << endl << endl;
  bool converged = false;
  for (int iter = 0; iter != MAX_ITER_FCI; ++iter) { 
    int start = ::clock();

    // form a sigma vector given cc
    form_sigma(cc, sigma, d, e, Jop);

// TODO TODO TODO -> this part is not multi-state 
// perhaps split set and get functions in davidson.h
    const vector<double> energies = davidson.compute(shared_ptr<Civec>(new Civec(*cc->data(0))),
                                                     shared_ptr<Civec>(new Civec(*sigma->data(0))));

    // get residual and new vectors
    vector<shared_ptr<Civec> > errvec = davidson.residual();

    // compute errors
    vector<double> errors;
    for (int i = 0; i != num_state; ++i) errors.push_back(errvec[i]->norm());

    if (*max_element(errors.begin(), errors.end()) < THREASH) {
      converged = true;
    } else { 
      // denominator scaling 
      for (int ist = 0; ist != num_state; ++ist) {
        const int size = cc->data(ist)->size();
        double* target_array = cc->data(ist)->first();
        double* source_array = errvec[ist]->first();
        double* denom_array = denom_->first();
        const double en = energies[ist];
        for (int i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / (en - denom_array[i]);
        }
  // TODO this will be changed to "add function in the future"
        davidson.orthog(cc->data(ist));
      }
    }

    // printing out
    int end = ::clock();
    for (int i = 0; i != num_state; ++i) {
      cout << indent << setw(5) << iter << setw(5) << i << setw(20) << fixed << setprecision(12) << energies[i]+nuclear << space3 
                                        << setw(17) << errors[i] << setw(15) << setprecision(2) << (end - start) * 1.0e-6 << endl; 
    }
    if (converged) break;
  }
  // main iteration ends here
}


void FCI::form_sigma(shared_ptr<Dvec> ccvec, shared_ptr<Dvec> sigmavec,
                     shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<MOFile> Jop) { // d and e are scratch area for D and E intermediates 
  const int la = d->lena();  
  const int lb = d->lenb();  
  const int ij = d->ij();
  sigmavec->zero();

  const int nstate = ccvec->ij();

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    {
      for (auto iter = phia_.begin();  iter != phia_.end(); ++iter) {
        const double hc = Jop->mo1e(get<2>(*iter)) * get<1>(*iter);
        daxpy_(&lb, &hc, cc->element_ptr(0, get<3>(*iter)), &unit, sigma->element_ptr(0, get<0>(*iter)), &unit); 
      }
    }
    // (task2) two electron contributions
    {
      // zero out intermediates.
      d->zero();

      // step (c)
      // (task2a-1) D(Phib, Phia, ij) += sign C(Psib, Phi'a)
      {
        for (auto iter = phia_.begin();  iter != phia_.end(); ++iter) {
          const double sign = static_cast<double>(get<1>(*iter));
          const int ip = get<2>(*iter);
          double* const target_array = d->data(ip)->element_ptr(0, get<3>(*iter));
          daxpy_(&lb, &sign, cc->element_ptr(0, get<0>(*iter)), &unit, target_array, &unit);
        }
      }

      // step (d)
      // (task2a-2) D(Phib, Phia, ij) += sign C(Psib', Phia)
      {
        for (int i = 0; i != la; ++i) {
          double* const source_array = cc->element_ptr(0, i); 
          for (int ip = 0; ip != ij; ++ip) {
            double* const target_array = d->data(ip)->element_ptr(0, i); 
            for (auto iter = phib_.begin(); iter != phib_.end(); ++iter) {
              if (get<2>(*iter) == ip) { // TODO dislike "if" here, but this seems the only way.. see (g) as well 
                const double sign = static_cast<double>(get<1>(*iter));
                target_array[get<3>(*iter)] += sign * source_array[get<0>(*iter)]; 
              }
            }
          } 
        }
      }

      // step (e)
      // (task2b) E(Phib, Phia, ij) = D(Psib, Phia, ij) (ij|kl)
      {
        const int lenab = la*lb;
        dgemm_("n", "n", &lenab, &ij, &ij, &half, d->first(), &lenab, Jop->mo2e_ptr(), &ij,
                                          &zero, e->first(), &lenab);
      }

      // step (f)
      // (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, ij)
      {
        for (auto iter = phia_.begin(); iter != phia_.end(); ++iter) {
          const double sign = static_cast<double>(get<1>(*iter)); 
          const int ip = get<2>(*iter);
          double* const target_array = sigma->element_ptr(0, get<0>(*iter));
          daxpy_(&lb, &sign, e->data(ip)->element_ptr(0, get<3>(*iter)), &unit, target_array, &unit);
        }
      }

      // step (g)
      // (task2c-2) sigma(Phib', Phia) += sign E(Psib, Phia, ij)
      {
        for (int i = 0; i != la; ++i) {
          double* const target_array = sigma->element_ptr(0, i);
          for (int ip = 0; ip != ij; ++ip) {
            double* const source_array = e->data(ip)->element_ptr(0, i);
            for (auto iter = phib_.begin(); iter != phib_.end(); ++iter) {
              if (get<2>(*iter) == ip) {
                const double sign = static_cast<double>(get<1>(*iter));
                target_array[get<0>(*iter)] += sign * source_array[get<3>(*iter)]; 
              }
            }
          }
        }
      }
    }

    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    {
      for (int i = 0; i != la; ++i) {
        double* const target_array = sigma->element_ptr(0, i);
        double* const source_array = cc->element_ptr(0, i);
        for (auto iter = phib_.begin();  iter != phib_.end(); ++iter) {
          const double hc = Jop->mo1e(get<2>(*iter)) * get<1>(*iter);
          target_array[get<0>(*iter)] += hc * source_array[get<3>(*iter)];
        }
      }
    }
  }
}


void FCI::const_denom(shared_ptr<MOFile> Jop) {

  vector<double> jop, kop;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*Jop->mo2e(j*norb_+j,i*norb_+i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*Jop->mo2e(i*norb_+j,i*norb_+j);
    }
  }
  shared_ptr<Civec> tmp(new Civec(stringb_.size(), stringa_.size()));
  denom_ = tmp;

  double* iter = denom_->first();
  for (auto ia = stringa_.begin(); ia != stringa_.end(); ++ia) {
    for (auto ib = stringb_.begin(); ib != stringb_.end(); ++ib, ++iter) {
      unsigned int iabit1 = *ia;
      unsigned int ibbit1 = *ib;
      for (int i = 0; i != norb_; ++i, (iabit1 >>= 1), (ibbit1 >>= 1)) {
        const unsigned int nia = (iabit1&1);
        const unsigned int nib = (ibbit1&1);
        *iter += Jop->mo1e(i*norb_+i) * (nia + nib);
        unsigned int iabit2 = *ia;
        unsigned int ibbit2 = *ib;
        for (int j = 0; j != norb_; ++j, (iabit2 >>= 1), (ibbit2 >>= 1)) {
          const unsigned int nja = (iabit2&1);
          const unsigned int njb = (ibbit2&1);
          const unsigned int addj = ((nia & njb) << 1) + (nia & nja) + (nib & njb);
          *iter += jop[j+norb_*i] * addj;
          const unsigned int addk = (nia & (1^nja)) + (nib & (1^njb));
          *iter += kop[j+norb_*i] * addk;
        }
      }
    }
  }
}


void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}
