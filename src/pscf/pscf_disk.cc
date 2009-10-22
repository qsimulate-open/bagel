//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <src/pscf/pscf_disk.h>
#include <src/rysint/eribatch.h>
#include <src/util/pcompfile.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<PGeometry> RefPGeometry;
typedef boost::shared_ptr<PFock> RefPFock;
typedef boost::shared_ptr<PTildeX> RefPTildeX;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<PMatrix1e> RefPMatrix1e;

using namespace std;
using namespace boost;

PSCF_DISK::PSCF_DISK(const RefPGeometry g) : PSCF(g) {

  store_ERI();

}

PSCF_DISK::~PSCF_DISK() {

}


void PSCF_DISK::store_ERI() {

  shared_ptr<PCompFile<ERIBatch> > tmpae(new PCompFile<ERIBatch>(geom_, 0.0, false, "ERI OBS"));
  ao_eri_ = tmpae;
  ao_eri_->store_integrals();
  ao_eri_->reopen_with_inout();

}

