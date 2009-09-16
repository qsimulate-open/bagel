//
// Author: Toru Shiozaki
// Date  : May 2009
//

// compute overlap and kinetic integrals using Obara-Saika recursion formula

#include <src/osint/overlapbatch.h>

using namespace std;

OverlapBatch::OverlapBatch(const vector<boost::shared_ptr<Shell> >& _basis) 
 : OSInt(_basis) {

}  


OverlapBatch::~OverlapBatch() {
}

