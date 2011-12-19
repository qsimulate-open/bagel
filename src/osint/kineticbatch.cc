//
// Author: Toru Shiozaki
// Date  : May 2009
//

// compute overlap and kinetic integrals using Obara-Saika recursion formula

#include <src/osint/kineticbatch.h>

using namespace std;

KineticBatch::KineticBatch(const vector<std::shared_ptr<Shell> >& _basis) 
 : OSInt(_basis) {

}  


KineticBatch::~KineticBatch() {
}

