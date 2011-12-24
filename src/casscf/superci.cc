//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/superci.h>
#include <iostream>
#include <src/fci/fci.h>

using namespace std;

SuperCI::SuperCI(const std::shared_ptr<Geometry> geom) : CASSCF(geom) {
  cout << "    * using the Super CI algorithm as noted in Roos (1980) IJQC" << endl << endl;

}

SuperCI::~SuperCI() {

}

void SuperCI::update_orbitals_() {

}

void SuperCI::update_civectors_() {

}

void SuperCI::compute() {
  // first compute reference.
  ref_->compute();

  // first perform CASCI
  streambuf* backup = cout.rdbuf();
  streambuf* null;
  cout.rdbuf(null);

  shared_ptr<FCI> fci(new FCI(geom_));

  cout.rdbuf(backup);


}
