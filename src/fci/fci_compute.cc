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

// TODO hardwired MAXIT.
#define MAX_ITER_FCI 1

using namespace std;

static const int unit = 1;
static const double one = 1.0;
static const double zero = 0.0; 
static const string indent = "  ";
static const string space3 = "   "; 

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

  // some constants
  const int la = stringa_.size();
  const int lb = stringb_.size();
  const int ij = norb_ * norb_; 

  // we need two vectors for intermediate quantities
  shared_ptr<Dvec> d(new Dvec(lb, la, ij));
  shared_ptr<Dvec> e(new Dvec(lb, la, ij));

  // Creating an initial CI vector
  shared_ptr<Civec> cc(new Civec(lb, la)); // B runs first
  shared_ptr<Civec> sigma(new Civec(lb, la));

  // TODO This is only if RHF is valid and therefore wrong in general cases.
  //      The right ways is to compute diagonal dinominators and select determinants
  //      that have small values, and then spin adapt.
  cc->element(0,0) = 1.0;

  // nuclear energy retrieved from geometry
  const double nuclear = geom_->nuclear_repulsion()

  // Davidson utility
  const int num_state = 1; // TODO should be read from the input
  DavidsonDiag<Civec> davidson(num_state, MAX_ITER_FCI);

  // main iteration starts here
  for (int iter = 0; iter != MAX_ITER_FCI; ++iter) { 
    int start = ::clock();

    double error = 0.0;

    form_sigma(cc, sigma, d, e, Jop);
cout << setprecision(12) << fixed;
cout << cc->norm() << endl;
cout << sigma->norm() << endl;
cout << cc->ddot(*sigma) << endl;
    const vector<double> energies = davidson.compute(cc, sigma);

    int end = ::clock();
    for (int i = 0; i != num_state; ++i) {
      cout << indent << setw(5) << iter << setw(5) << i << setw(20) << fixed << setprecision(12) << energies[i] << space3 
                                        << setw(17) << error << setw(15) << setprecision(2) << (end - start) * 1.0e-6 << endl; 
    }
  }
  // main iteration ends here
}


void FCI::form_sigma(shared_ptr<Civec> cc, shared_ptr<Civec> sigma,
                     shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<MOFile> Jop) { // d and e are scratch area for D and E intermediates 
  const int la = d->lena();  
  const int lb = d->lenb();  
  const int ij = d->ij();
  sigma->zero();

  // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
  {
    for (auto iter = phia_.begin();  iter != phia_.end(); ++iter) {
      const double hc = Jop->mo1e(get<2>(*iter)) * get<1>(*iter);
      daxpy_(&lb, &hc, cc->element_ptr(0, get<3>(*iter)), &unit, sigma->element_ptr(0, get<0>(*iter)), &unit); 
    }
  }
#if 1
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
      dgemm_("n", "n", &lenab, &ij, &ij, &one, d->first(), &lenab, Jop->mo2e_ptr(), &ij,
                                        &zero, e->first(), &lenab);
    }

    // step (f)
    // (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, ij)
    for (auto iter = phia_.begin(); iter != phia_.end(); ++iter) {
      const double sign = static_cast<double>(get<1>(*iter)); 
      const int ip = get<2>(*iter);
      double* const target_array = sigma->element_ptr(0, get<0>(*iter));
      daxpy_(&lb, &sign, e->data(ip)->element_ptr(0, get<3>(*iter)), &unit, target_array, &unit);
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
#endif

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

void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}
