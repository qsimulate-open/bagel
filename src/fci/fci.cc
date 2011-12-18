//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <iostream>
#include <stdexcept>
#include <src/fci/mofile.h>
#include <src/fci/fci.h>
#include <src/fci/comb.h>
#include <boost/algorithm/combination.hpp>

using namespace boost;
using namespace std;

static const Comb comb;

FCI::FCI(const shared_ptr<Geometry> geom)
 : ref_(new SCF(geom)), geom_(geom), nelea_(geom_->nocc()/2), neleb_(geom_->nocc()/2),
   ncore_(0), norb_(geom_->nbasis()) {
  // ^- TODO somehow we need the input interface to number of eletrons in alpha and beta!!!

  const_lexical_mapping_();
  const_string_lists_();
  const_phis_();
}

FCI::~FCI() {

}


void FCI::const_phis_() {
  typedef vector<tuple<unsigned int, int, unsigned int, unsigned int> >::iterator Titer;

  phia_.resize(stringa_.size()*norb_*norb_);
  Titer piter = phia_.begin();

  for (vector<unsigned int>::const_iterator iter = stringa_.begin(); iter != stringa_.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { // annihilation
      const unsigned int ibit = (1 << i);
      if (ibit & *iter) {
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & *iter)) {
            const unsigned int mbit = jbit^*iter;
            *piter = make_tuple(lexicala(mbit), sign(mbit, i, j), j, i); 
            ++piter;
          }
        }
      }
    }
  }

}

void FCI::const_string_lists_() {
  vector<int> data(norb_);
  for (int i=0; i!=norb_; ++i)  data[i] = i;

  const int lengtha = comb.c(norb_, nelea_);
  const int lengthb = comb.c(norb_, neleb_);
  stringa_.resize(lengtha); 
  stringb_.resize(lengthb);
  fill(stringa_.begin(), stringa_.end(), 0);
  fill(stringb_.begin(), stringb_.end(), 0);

  vector<unsigned int>::iterator sa = stringa_.begin(); 
  do {
    for (int i=0; i!=nelea_; ++i) *sa += (1 << data[i]);
    ++sa;
  } while (next_combination(data.begin(), data.begin()+nelea_, data.end()));

  sa = stringb_.begin(); 
  do {
    for (int i=0; i!=neleb_; ++i) *sa += (1 << data[i]);
    ++sa;
  } while (next_combination(data.begin(), data.begin()+neleb_, data.end()));

#if 0
  for (vector<unsigned int>::const_iterator i = stringa_.begin(); i != stringa_.end(); ++i) {
    cout << lexicala(*i) << endl;
  }
  for (vector<unsigned int>::const_iterator i = stringb_.begin(); i != stringb_.end(); ++i) {
    cout << lexicalb(*i) << endl;
  }
#endif

}

void FCI::const_lexical_mapping_() {
  // combination numbers up to 31 orbitals (fci/comb.h)
  zkla_.resize(nelea_ * norb_); 
  zklb_.resize(neleb_ * norb_);
  fill(zkla_.begin(), zkla_.end(), 0llu);
  fill(zklb_.begin(), zklb_.end(), 0llu);

  // this part is 1 offset due to the convention of Knowles & Handy's paper.
  for (int k = 1; k < nelea_; ++k) {
    for (int l = k; l <= norb_-nelea_+k; ++l) {
      for (int m = norb_-l+1; m <= norb_-k; ++m) {
        zkla(k-1, l-1) += comb.c(m, nelea_-k) - comb.c(m-1, nelea_-k-1); 
      }
    }
  }
  for (int l = nelea_; l <= norb_; ++l) zkla(nelea_-1, l-1) = l - nelea_; 

  if (nelea_ == neleb_) {
    copy(zkla_.begin(), zkla_.end(), zklb_.begin());
  } else {
    for (int k = 1; k <= neleb_; ++k) {
      for (int l = k; l <= norb_-neleb_+k; ++l) {
        zklb(k-1, l-1) = 0; 
        for (int m = norb_-l+1; m <= norb_-k; ++m) {
          zklb(k-1, l-1) += comb.c(m, neleb_-k) - comb.c(m-1, neleb_-k-1); 
        }
      }
    }
    for (int l = neleb_; l <= norb_; ++l) zklb(neleb_-1, l-1) = l - neleb_; 
  }
}

void FCI::compute() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("FCI: C1 only at the moment."); 

  cout << "  ---------------------------" << endl;
  cout << "  |     FCI calculation     |" << endl;
  cout << "  ---------------------------" << endl << endl;

  // first obtain reference function
  ref_->compute();
  shared_ptr<Coeff> cmo = ref_->coeff(); 

  // iiii file to be created (MO transformation).
  // now Jop->mo1e() and Jop->mo2e() contains one and two body part of Hamiltonian
  shared_ptr<MOFile> Jop(new MOFile(geom_, cmo));

  // TODO right now full basis is used. 
  Jop->create_Jiiii(ncore_, ncore_+norb_);

  // note h'_kl = h_kl - 0.5 (k*|*l)
  // TODO

  // (alpha beta) two excitations

}
