//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <iostream>
#include <iomanip>
#include <src/fci/fci.h>
#include <src/fci/comb.h>
#include <boost/algorithm/combination.hpp>

// TODO hardwired
#include <src/fci/macros.h>

using namespace boost;
using namespace std;

static const Comb comb;

FCI::FCI(const std::shared_ptr<Geometry> geom)
 : ref_(new SCF(geom)), geom_(geom), nelea_(geom_->nocc()/2 - geom->nfrc()/2), neleb_(geom_->nocc()/2 - geom->nfrc()/2),
   ncore_(geom->nfrc()/2), norb_(geom_->nbasis() - geom->nfrc()/2), num_state_(NUM_STATE) {
  // ^- TODO somehow we need the input interface to number of eletrons in alpha and beta!!!
  ref_->compute();
  common_init();
}

FCI::FCI(const std::shared_ptr<Geometry> geom, shared_ptr<SCF> r)
 : ref_(r), geom_(geom), nelea_(geom_->nocc()/2 - geom->nfrc()/2), neleb_(geom_->nocc()/2 - geom->nfrc()/2),
   ncore_(geom->nfrc()/2), norb_(geom_->nbasis() - geom->nfrc()/2), num_state_(NUM_STATE) {
  // ^- TODO somehow we need the input interface to number of eletrons in alpha and beta!!!
  common_init();
}


void FCI::common_init() {
  print_header();

  cout << "  Performs exactly the same way as Knowles & Handy 1984 CPL" << endl;
  cout << endl;
  cout << "  o lexical mappings" << endl;
  const_lexical_mapping_();
  cout << "  o alpha-beta strings" << endl;
  const_string_lists_();
  cout << "      length: " << setw(13) << stringa_.size() + stringb_.size() << endl;
  cout << "  o single displacement lists (alpha)" << endl;
  const_phis_<0>(stringa_, phia_);
  cout << "      length: " << setw(13) << phia_.size()*phia_.front().size() << endl;
  cout << "  o single displacement lists (beta)" << endl;
  const_phis_<1>(stringb_, phib_);
  cout << "      length: " << setw(13) << phib_.size()*phib_.front().size() << endl;
  cout << endl;

}

FCI::~FCI() {

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

  auto sa = stringa_.begin(); 
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
  for (auto i = stringa_.begin(); i != stringa_.end(); ++i) {
    cout << lexical<0>(*i) << endl;
  }
  for (auto i = stringb_.begin(); i != stringb_.end(); ++i) {
    cout << lexical<1>(*i) << endl;
  }
#endif

}

void FCI::const_lexical_mapping_() {
  // combination numbers up to 31 orbitals (fci/comb.h)
  zkl_.resize(nelea_ * norb_ + neleb_ * norb_); 
  fill(zkl_.begin(), zkl_.end(), 0u);

  // this part is 1 offset due to the convention of Knowles & Handy's paper.
  // Just a blind copy from the paper without understanding much, but the code below works. 
  for (int k = 1; k < nelea_; ++k) {
    for (int l = k; l <= norb_-nelea_+k; ++l) {
      for (int m = norb_-l+1; m <= norb_-k; ++m) {
        zkl(k-1, l-1, Alpha) += comb.c(m, nelea_-k) - comb.c(m-1, nelea_-k-1); 
      }
    }
  }
  for (int l = nelea_; l <= norb_; ++l) zkl(nelea_-1, l-1, Alpha) = l - nelea_; 

  if (nelea_ == neleb_) {
    copy(zkl_.begin(), zkl_.begin() + nelea_*norb_, zkl_.begin() + nelea_*norb_);
  } else {
    for (int k = 1; k <= neleb_; ++k)
      for (int l = k; l <= norb_-neleb_+k; ++l)
        for (int m = norb_-l+1; m <= norb_-k; ++m)
          zkl(k-1, l-1, Beta) += comb.c(m, neleb_-k) - comb.c(m-1, neleb_-k-1); 
    for (int l = neleb_; l <= norb_; ++l) zkl(neleb_-1, l-1, Beta) = l - neleb_; 
  }
}


void FCI::print_civectors(vector<shared_ptr<Civec> > vec, const double thr) const {
  int j = 0;
  for (auto iter = vec.begin(); iter != vec.end(); ++iter, ++j) {
    cout << endl;
    cout << "     * ci vector, state " << setw(3) << j << endl; 
    const double* i = (*iter)->first();
    multimap<double, tuple<double, int, int> > tmp;
    for (auto ia = stringa_.begin(); ia != stringa_.end(); ++ia) {
      for (auto ib = stringb_.begin(); ib != stringb_.end(); ++ib, ++i) {
        if (abs(*i) > thr) {
          tmp.insert(make_pair(-abs(*i), make_tuple(*i, *ia, *ib))); 
        }
      }
    }
    for (auto iter = tmp.begin(); iter != tmp.end(); ++iter) {
      cout << "       " << print_bit(get<1>(iter->second), get<2>(iter->second))
           << "  " << setprecision(10) << setw(15) << get<0>(iter->second) << endl; 
    }
  }
}


void FCI::print_timing_(const string label, int& time, std::vector<pair<string, double> >& timing) const {
  timing.push_back(make_pair(label, (::clock()-time)/static_cast<double>(CLOCKS_PER_SEC)));
  time = ::clock();
}
