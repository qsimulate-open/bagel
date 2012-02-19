//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/storage.h>
#include <src/smith/tensor.h>
#include <iostream>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/smith/moint.h>
#include <src/wfn/reference.h>

using namespace SMITH;
using namespace std;

void a(shared_ptr<Reference> r){
  const int max = 7;
  IndexRange closed(10, max);
  IndexRange acc(7, max);
  IndexRange virt(20, max);

#if 0
  closed.print();
  acc.print();
  virt.print();
#endif

  vector<IndexRange> o;
  o.push_back(closed);
  o.push_back(virt);
  o.push_back(closed);
  o.push_back(virt);

  MOInt<Storage_Incore> a(r, o);
  shared_ptr<Tensor<Storage_Incore> > tensor = a.data();

  // debug implementation of MP2 here.
  shared_ptr<Fock<1> > fock0(new Fock<1>(r->geom(), r->hcore()));
  shared_ptr<Matrix1e> den(new Matrix1e(r->coeff()->form_density_rhf()));
  shared_ptr<Fock<1> > fock1(new Fock<1>(r->geom(), fock0, den, r->schwarz()));
  Matrix1e f = *r->coeff() % *fock1 * *r->coeff();

  vector<double> eig;
  const int nb = r->nclosed() + r->nact() + r->nvirt();
  for (int i = 0; i != nb; ++i) { eig.push_back(f.element(i,i)); }

  for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
    for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
      for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
        for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
          vector<size_t> h;
          h.push_back(i0->key());
          h.push_back(i1->key());
          h.push_back(i2->key());
          h.push_back(i3->key());
          unique_ptr<double[]> d = tensor->get_block(h);
        }
      }
    }
  }

}
