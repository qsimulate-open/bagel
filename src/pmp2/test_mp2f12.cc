//
// Author : Toru Shiozaki
// Date   : April 27
//

#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <src/pscf/pgeometry.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/pscf.h>
#include <src/pscf/pscf_disk.h>
#include <src/pmp2/pmp2.h>

#include <src/stackmem.h>

#include <src/util/input.h>

using namespace std;
extern  StackMem* stack;

void test_mp2f12() {
  StackMem* a = new StackMem(static_cast<size_t>(100000000LU));
  stack = a;

  shared_ptr<PGeometry> pgeom(new PGeometry("oldinp", 0));
  shared_ptr<PSCF_DISK> pscf(new PSCF_DISK(pgeom));
  pscf->compute();

  PMP2 pmp2(pscf->geom(), pscf->coeff(), pscf->eig(), pscf->ao_eri(), false);
  pmp2.compute();


};

