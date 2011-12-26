//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/casscf.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

using namespace std;

CASSCF::CASSCF(multimap<string, string> idat, const shared_ptr<Geometry> geom) :
  idata_(idat), geom_(geom), ref_(new SCF(idat, geom)) {

  ref_->compute();
  common_init();

};

CASSCF::CASSCF(multimap<string, string> idat, const shared_ptr<Geometry> geom, const shared_ptr<SCF> scf)
  : idata_(idat), geom_(geom), ref_(scf) {

  common_init();

};


void CASSCF::common_init() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment."); 
  print_header();

  {// get maxiter from the input
    max_iter_ = 100;
    auto iter = idata_.find("maxiter");
    if (iter != idata_.end()) max_iter_ = boost::lexical_cast<int>(iter->second);
  }
  {// get nstate from the input
    nstate_ = 1;
    auto iter = idata_.find("nstate");
    if (iter != idata_.end()) nstate_ = boost::lexical_cast<int>(iter->second);
  }

  // nocc from the input. If not present, full valence active space is generated.
  {
    auto iter = idata_.find("nocc");
    if (iter != idata_.end()) {
      nocc_ = boost::lexical_cast<int>(iter->second);
      if (nocc_ <= 0) throw runtime_error("It appears that nocc <= 0. Check nocc value.");
    } else {
      cout << "    * full valence occupied space generated for nocc." << endl;
      nocc_ = geom_->num_count_full_valence_nocc(); 
    }
  }
  // nclosed from the input. If not present, full core space is generated.
  {
    auto iter = idata_.find("nclosed");
    if (iter != idata_.end()) {
      nclosed_ = boost::lexical_cast<int>(iter->second);
    } else {
      nclosed_ = geom_->num_count_ncore() / 2; 
      cout << "    * full core space generated for nclosed." << endl;
    }
  }
  nact_ = nocc_ - nclosed_;
  if (nact_ <= 0) throw runtime_error("It appears that nact <= 0. Check the nocc and nclosed values");
  nbasis_ = geom_->nbasis();
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;


  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and total orbitals.
  mute_stdcout();
  {
    multimap<string, string> fullci_data = idata_;
    for (auto iter = fullci_data.equal_range("ncore").first; iter != fullci_data.equal_range("ncore").second; ++iter)
      fullci_data.erase(iter);
    fullci_data.insert(make_pair("ncore", boost::lexical_cast<string>(nclosed_)));
    for (auto iter = fullci_data.equal_range("norb").first; iter != fullci_data.equal_range("norb").second; ++iter)
      fullci_data.erase(iter);
    fullci_data.insert(make_pair("norb", boost::lexical_cast<string>(nact_)));
    shared_ptr<FCI> fci_tmp(new FCI(fullci_data, geom_, ref_));
    fci_ = fci_tmp;
  }
//fci_->compute();
  resume_stdcout();

};


CASSCF::~CASSCF() {

};


void CASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}

static streambuf* backup_stream_;

void CASSCF::mute_stdcout() {
  backup_stream_ = cout.rdbuf();
  streambuf* null;
  cout.rdbuf(null);
}


void CASSCF::resume_stdcout() {
  cout.rdbuf(backup_stream_);
}

