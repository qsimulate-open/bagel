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

static int address(int i, int j) { assert(i <= j); return i+((j*(j+1))>>1); };

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
  const int ij = norb_ * (norb_+1) /2;

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
          target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
        }
  // TODO this will be changed to "add function in the future"
        davidson.orthog(cc->data(ist));
      }
    }

    // printing out
    int end = ::clock();
    for (int i = 0; i != num_state; ++i) {
      cout << indent << setw(5) << iter << setw(5) << i << setw(20) << fixed << setprecision(12) << energies[i]+nuclear << space3 
                                        << setw(17) << errors[i] << setw(15) << setprecision(2)
                                        << (end - start)/static_cast<double>(CLOCKS_PER_SEC) << endl; 
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
//#define TIMING

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
#ifdef TIMING
    vector<pair<string, double> > timing;
    int start = ::clock();
#endif
    {
      for (int ip = 0; ip != ij; ++ip) {
        const double h = Jop->mo1e(ip);
        for (auto iter = phia_[ip].begin();  iter != phia_[ip].end(); ++iter) {
          const double hc = h * get<1>(*iter);
          daxpy_(&lb, &hc, cc->element_ptr(0, get<2>(*iter)), &unit, sigma->element_ptr(0, get<0>(*iter)), &unit); 
        }
      }
    }
#ifdef TIMING
    { const string task("task1"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif
    // (task2) two electron contributions
    {
      // zero out intermediates.
      d->zero();

      // step (c)
      // (task2a-1) D(Phib, Phia, ij) += sign C(Psib, Phi'a)
      {
        double* const source_base = cc->first();
        for (int ip = 0; ip != ij; ++ip) {
          double* const target_base = d->data(ip)->first();
          for (auto iter = phia_[ip].begin();  iter != phia_[ip].end(); ++iter) {
            const double sign = static_cast<double>(get<1>(*iter));
            double* const target_array = target_base + get<2>(*iter)*lb;
            daxpy_(&lb, &sign, source_base + get<0>(*iter)*lb, &unit, target_array, &unit);
          }
        }
      }
#ifdef TIMING
    { const string task("task2a-1"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif

      // step (d)
      // (task2a-2) D(Phib, Phia, ij) += sign C(Psib', Phia)
      {
        for (int i = 0; i < la; i+=4) {
          if (i+3 < la) {
            double* const source_array0 = cc->element_ptr(0, i); 
            double* const source_array1 = cc->element_ptr(0, i+1); 
            double* const source_array2 = cc->element_ptr(0, i+2); 
            double* const source_array3 = cc->element_ptr(0, i+3); 
            for (int ip = 0; ip != ij; ++ip) {
              double* const target_array0 = d->data(ip)->element_ptr(0, i); 
              double* const target_array1 = d->data(ip)->element_ptr(0, i+1); 
              double* const target_array2 = d->data(ip)->element_ptr(0, i+2); 
              double* const target_array3 = d->data(ip)->element_ptr(0, i+3); 
              for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
                const double sign = static_cast<double>(get<1>(*iter));
                const int ia = get<2>(*iter);
                const int ib = get<0>(*iter);
                target_array0[ia] += sign * source_array0[ib]; 
                target_array1[ia] += sign * source_array1[ib]; 
                target_array2[ia] += sign * source_array2[ib]; 
                target_array3[ia] += sign * source_array3[ib]; 
              }
            } 
          } else {
            for (int j = i; j != la; ++j) {
              double* const source_array0 = cc->element_ptr(0, j);
              for (int ip = 0; ip != ij; ++ip) {
                double* const target_array0 = d->data(ip)->element_ptr(0, j);
                for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
                  const double sign = static_cast<double>(get<1>(*iter));
                  target_array0[get<2>(*iter)] += sign * source_array0[get<0>(*iter)]; 
                }
              }
            }
          }
        }
      }
#ifdef TIMING
    { const string task("task2a-2"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif

      // step (e)
      // (task2b) E(Phib, Phia, ij) = D(Psib, Phia, ij) (ij|kl)
      {
        const int lenab = la*lb;
        dgemm_("n", "n", &lenab, &ij, &ij, &half, d->first(), &lenab, Jop->mo2e_ptr(), &ij,
                                           &zero, e->first(), &lenab);
      }
#ifdef TIMING
    { const string task("task2b"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif

      // step (f)
      // (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, ij)
      {
        for (int ip = 0; ip != ij; ++ip) { 
          double* const source_base = e->data(ip)->first();
          for (auto iter = phia_[ip].begin(); iter != phia_[ip].end(); ++iter) {
            const double sign = static_cast<double>(get<1>(*iter)); 
            double* const target_array = sigma->element_ptr(0, get<0>(*iter));
            daxpy_(&lb, &sign, source_base + lb*get<2>(*iter), &unit, target_array, &unit);
          }
        }
      }
#ifdef TIMING
    { const string task("task2c-1"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif

      // step (g)
      // (task2c-2) sigma(Phib', Phia) += sign E(Psib, Phia, ij)
      {
        for (int i = 0; i < la; i+=8) {
        if (i+7 < la) {
          double* const target_array0 = sigma->element_ptr(0, i);
          double* const target_array1 = sigma->element_ptr(0, i+1);
          double* const target_array2 = sigma->element_ptr(0, i+2);
          double* const target_array3 = sigma->element_ptr(0, i+3);
          double* const target_array4 = sigma->element_ptr(0, i+4);
          double* const target_array5 = sigma->element_ptr(0, i+5);
          double* const target_array6 = sigma->element_ptr(0, i+6);
          double* const target_array7 = sigma->element_ptr(0, i+7);
          for (int ip = 0; ip != ij; ++ip) {
            double* const source_array0 = e->data(ip)->element_ptr(0, i);
            double* const source_array1 = e->data(ip)->element_ptr(0, i+1);
            double* const source_array2 = e->data(ip)->element_ptr(0, i+2);
            double* const source_array3 = e->data(ip)->element_ptr(0, i+3);
            double* const source_array4 = e->data(ip)->element_ptr(0, i+4);
            double* const source_array5 = e->data(ip)->element_ptr(0, i+5);
            double* const source_array6 = e->data(ip)->element_ptr(0, i+6);
            double* const source_array7 = e->data(ip)->element_ptr(0, i+7);
            for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
              const double sign = static_cast<double>(get<1>(*iter));
              target_array0[get<0>(*iter)] += sign * source_array0[get<2>(*iter)]; 
              target_array1[get<0>(*iter)] += sign * source_array1[get<2>(*iter)]; 
              target_array2[get<0>(*iter)] += sign * source_array2[get<2>(*iter)]; 
              target_array3[get<0>(*iter)] += sign * source_array3[get<2>(*iter)]; 
              target_array4[get<0>(*iter)] += sign * source_array4[get<2>(*iter)]; 
              target_array5[get<0>(*iter)] += sign * source_array5[get<2>(*iter)]; 
              target_array6[get<0>(*iter)] += sign * source_array6[get<2>(*iter)]; 
              target_array7[get<0>(*iter)] += sign * source_array7[get<2>(*iter)]; 
            }
          }
        } else {
          for (int j = i; j != la; ++j) {
            double* const target_array0 = sigma->element_ptr(0, j);
            for (int ip = 0; ip != ij; ++ip) {
              double* const source_array0 = e->data(ip)->element_ptr(0, j);
              for (auto iter = phib_[ip].begin(); iter != phib_[ip].end(); ++iter) {
                const double sign = static_cast<double>(get<1>(*iter));
                target_array0[get<0>(*iter)] += sign * source_array0[get<2>(*iter)]; 
              }
            }
          }
        }
        }
      }
    }
#ifdef TIMING
    { const string task("task2c-2"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
#endif

    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    {
      for (int i = 0; i < la; i+=8) {
        if (i+7 < la) {
          double* const target_array0 = sigma->element_ptr(0, i);
          double* const target_array1 = sigma->element_ptr(0, i+1);
          double* const target_array2 = sigma->element_ptr(0, i+2);
          double* const target_array3 = sigma->element_ptr(0, i+3);
          double* const target_array4 = sigma->element_ptr(0, i+4);
          double* const target_array5 = sigma->element_ptr(0, i+5);
          double* const target_array6 = sigma->element_ptr(0, i+6);
          double* const target_array7 = sigma->element_ptr(0, i+7);
          double* const source_array0 = cc->element_ptr(0, i);
          double* const source_array1 = cc->element_ptr(0, i+1);
          double* const source_array2 = cc->element_ptr(0, i+2);
          double* const source_array3 = cc->element_ptr(0, i+3);
          double* const source_array4 = cc->element_ptr(0, i+4);
          double* const source_array5 = cc->element_ptr(0, i+5);
          double* const source_array6 = cc->element_ptr(0, i+6);
          double* const source_array7 = cc->element_ptr(0, i+7);
          for (int ip = 0; ip != ij; ++ip) {
            const double h = Jop->mo1e(ip);
            for (auto iter = phib_[ip].begin();  iter != phib_[ip].end(); ++iter) {
              const double hc = h * get<1>(*iter);
              target_array0[get<0>(*iter)] += hc * source_array0[get<2>(*iter)];
              target_array1[get<0>(*iter)] += hc * source_array1[get<2>(*iter)];
              target_array2[get<0>(*iter)] += hc * source_array2[get<2>(*iter)];
              target_array3[get<0>(*iter)] += hc * source_array3[get<2>(*iter)];
              target_array4[get<0>(*iter)] += hc * source_array4[get<2>(*iter)];
              target_array5[get<0>(*iter)] += hc * source_array5[get<2>(*iter)];
              target_array6[get<0>(*iter)] += hc * source_array6[get<2>(*iter)];
              target_array7[get<0>(*iter)] += hc * source_array7[get<2>(*iter)];
            }
          }
        } else {
          for (int j = i; j != la; ++j) {  
            double* const target_array0 = sigma->element_ptr(0, j);
            double* const source_array0 = cc->element_ptr(0, j);
            for (int ip = 0; ip != ij; ++ip) {
              const double h = Jop->mo1e(ip);
              for (auto iter = phib_[ip].begin();  iter != phib_[ip].end(); ++iter) {
                const double hc = h * get<1>(*iter);
                target_array0[get<0>(*iter)] += hc * source_array0[get<2>(*iter)];
              }
            }
          }
        }
      }
    }
#ifdef TIMING
    { const string task("task3"); timing.push_back(make_pair(task, (::clock()-start)/static_cast<double>(CLOCKS_PER_SEC))); start = ::clock(); }
    cout << "     timing info" << endl;
    for (auto iter = timing.begin(); iter != timing.end(); ++iter) {
      cout << "    " << setw(10) << iter->first << setw(10) << setprecision(2) << iter ->second << endl;
    }
#endif
  }
}


//
// averaged diagonal elements as defined in Knowles & Handy (1989) Compt. Phys. Comm. 
//
void FCI::const_denom(shared_ptr<MOFile> Jop) {

  vector<double> jop, kop, fk;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  fk.resize(norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*Jop->mo2e(j, j, i, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*Jop->mo2e(j, i, j, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    fk[i] = 0.0;
    for (int j = 0; j != norb_; ++j) {
      fk[i] += kop[i*norb_+j];
    }
  }
  shared_ptr<Civec> tmp(new Civec(stringb_.size(), stringa_.size()));
  denom_ = tmp;
  const int nspin = numofbits(stringa_.front()) - numofbits(stringb_.front());
  const int nspin2 = nspin*nspin;

  double* iter = denom_->first();
  for (auto ia = stringa_.begin(); ia != stringa_.end(); ++ia) {
    for (auto ib = stringb_.begin(); ib != stringb_.end(); ++ib, ++iter) {
      const int nopen = numofbits(iabit1^ibbit1);
      const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
      *iter = 0.0;
      unsigned int iabit1 = *ia;
      unsigned int ibbit1 = *ib;
      for (int i = 0; i != norb_; ++i, (iabit1 >>= 1), (ibbit1 >>= 1)) {
        const int nia = (iabit1&1);
        const int nib = (ibbit1&1);
        const int niab = nia + nib;
        const int Ni = (nia ^ nib);
        unsigned int iabit2 = *ia;
        unsigned int ibbit2 = *ib;
        for (int j = 0; j != i; ++j, (iabit2 >>= 1), (ibbit2 >>= 1)) {
          const int nja = (iabit2&1);
          const int njb = (ibbit2&1);
          const int Nj = (nja ^ njb);
          const int addj = niab * (nja + njb); 
          *iter += jop[j+norb_*i] * 2.0 * addj - kop[j+norb_*i] * (F*Ni*Nj + addj);
        }
        *iter += (Jop->mo1e(i,i) + fk[i]) * niab - kop[i+norb_*i] * 0.5 * (Ni - niab*niab);
      }
    }
  }
}


void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}
