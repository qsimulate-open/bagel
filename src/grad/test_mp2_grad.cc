//
// Author : Toru Shiozaki
// Date   : May 2012
//

// testing MP2 gradients

#include <src/wfn/reference.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace std;

// plain implementation of MP2 gradient to be replaced

void test_mp2_grad(shared_ptr<Reference> ref) {


  cout << "   testing MP2 grad implementation" << endl;

  assert(ref->nact() == 0);
  const size_t nocc = ref->nocc();
  const size_t nvirt = ref->nvirt();

  shared_ptr<const Geometry> geom = ref->geom();
  shared_ptr<const DensityFit> df = geom->df();

  unique_ptr<double[]> kijrs(new double[nocc*nvirt*nocc*nvirt]);


};
