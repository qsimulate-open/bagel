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
  idata_(idat), geom_(geom) {

  SCF scf(idat, geom);
  scf.compute();
  ref_ = scf.conv_to_ref();

  common_init();

};

CASSCF::CASSCF(multimap<string, string> idat, const shared_ptr<Geometry> geom, const shared_ptr<Reference> scf)
  : idata_(idat), geom_(geom), ref_(scf) {

  common_init();

};


void CASSCF::common_init() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment.");
  print_header();

  // get maxiter from the input
  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  // get maxiter from the input
  max_micro_iter_ = read_input<int>(idata_, "maxiter_micro", 100);
  // get nstate from the input
  nstate_ = read_input<int>(idata_, "nstate", 1);
  // get thresh (for macro iteration) from the input
  thresh_ = read_input<double>(idata_, "thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = read_input<double>(idata_, "thresh_micro", thresh_);

  // nocc from the input. If not present, full valence active space is generated.
  nocc_ = read_input<int>(idata_, "nocc", 0);
  if (nocc_ < 0) {
    throw runtime_error("It appears that nocc < 0. Check nocc value.");
  } else if (nocc_ == 0) {
    cout << "    * full valence occupied space generated for nocc." << endl;
    nocc_ = geom_->num_count_full_valence_nocc();
  }

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = read_input<int>(idata_, "nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore() / 2;
  }

  nact_ = nocc_ - nclosed_;
  if (nact_ <= 0) throw runtime_error("It appears that nact <= 0. Check the nocc and nclosed values");
#if 0
  nbasis_ = geom_->nbasis();
#else
  nbasis_ = ref_->coeff()->mdim();
#endif
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;


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
static ofstream* ofs_;

void CASSCF::mute_stdcout() {
  ofstream* ofs(new ofstream("casscf.log",(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
}


void CASSCF::resume_stdcout() {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


shared_ptr<Matrix1e> CASSCF::ao_rdm1(shared_ptr<RDM<1> > rdm1, const bool inactive_only) const {
  // first make 1RDM in MO
  shared_ptr<Matrix1e> mo_rdm1(new Matrix1e(geom_));
  for (int i = 0; i != nclosed_; ++i) mo_rdm1->element(i,i) = 1.0;
  if (!inactive_only) {
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nact_; ++j) {
        mo_rdm1->element(nclosed_+j, nclosed_+i) = rdm1->element(j,i)*0.5; // note the difference in notation between SCF and FCI...
      }
    }
  }
  // transform into AO basis
  const shared_ptr<Coeff> coeff = ref_->coeff();
  shared_ptr<Matrix1e> ao_rdm1(new Matrix1e(*coeff * *mo_rdm1 ^ *coeff)); // TODO make sure when overlap is truncated
  return ao_rdm1;
}

