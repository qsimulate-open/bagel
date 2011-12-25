//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/casscf.h>
#include <stdexcept>
#include <iostream>

using namespace std;

CASSCF::CASSCF(const std::shared_ptr<Geometry> geom) :
  geom_(geom), ref_(new SCF(geom)) {

  ref_->compute();

  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment."); 
  print_header();

  // CASSCF methods should have FCI member.
  shared_ptr<Geometry> local_geom(new Geometry(*geom_));
  // TODO input should be parsed.
  // TODO currently they are hard wired.
#if 0
  shared_ptr<FCI> fci_tmp(new FCI(local_geom, ref_));
  fci_ = fci_tmp;
fci_->compute();
#endif

};


CASSCF::~CASSCF() {

};


void CASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}
