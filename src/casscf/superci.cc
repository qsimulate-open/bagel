//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <src/casscf/superci.h>
#include <iostream>
#include <src/fci/fci.h>

using namespace std;

SuperCI::SuperCI(const multimap<string, string> idat, const shared_ptr<Geometry> geom)
 : CASSCF(idat, geom) {
  common_init();
}

SuperCI::SuperCI(const multimap<string, string> idat, const shared_ptr<Geometry> geom, const shared_ptr<SCF> scf)
 : CASSCF(idat, geom, scf) {
  common_init();
}

void SuperCI::common_init() {
  cout << "    * using the Super CI algorithm as noted in Roos (1980) IJQC" << endl << endl;
}

SuperCI::~SuperCI() {

}

void SuperCI::update_orbitals_() {

}

void SuperCI::update_civectors_() {

}

void SuperCI::compute() {

  // first perform CASCI


#if 0
#endif


}
