//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/hcore.h>
#include <src/osint/kineticbatch.h>
#include <src/rysint/naibatch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;

typedef boost::shared_ptr<Geometry> RefGeometry;
typedef boost::shared_ptr<Shell> RefShell;

Hcore::Hcore(const RefGeometry geom) : Matrix1e(geom) {

  init();

}


Hcore::~Hcore() {

}


void Hcore::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int nbasis) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis(); 
  const int dimb0 = input[1]->nbasis(); 
  KineticBatch kinetic(input);
  kinetic.compute();
  const double* kdata = kinetic.data();
  NAIBatch nai(input, geom_);
  nai.compute();
  const double* ndata = nai.data();

  int cnt = 0;
  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
      data_[i * nbasis + j] = kdata[cnt] + ndata[cnt];
    }
  }
}


