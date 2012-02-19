//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/prim_op.h>
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
  const int max = 10;
  IndexRange closed(r->nclosed(), max);
//IndexRange acc(r->nact(), max);
  IndexRange virt(r->nvirt(), max);

#if 0
  closed.print();
//acc.print();
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
  const int nocc = r->nclosed() + r->nact();
  const int nb = r->nclosed() + r->nact() + r->nvirt();
  for (int i = 0; i != nb; ++i) { eig.push_back(f.element(i,i)); }

  double en = 0.0;
  for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
    for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
      for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
        for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
          vector<size_t> h,g;
          h.push_back(i0->key()); h.push_back(i1->key()); h.push_back(i2->key()); h.push_back(i3->key());
          g.push_back(i0->key()); g.push_back(i3->key()); g.push_back(i2->key()); g.push_back(i1->key());
          const size_t size = tensor->get_size(h);
          assert(size == tensor->get_size(g));

          const unique_ptr<double[]> d = tensor->get_block(h);
          const unique_ptr<double[]> e = tensor->get_block(g);
          unique_ptr<double[]> buf(new double[size]);

          sort_indices4(e, buf, i1->size(), i2->size(), i3->size(), i0->size(), 2, 1, 0, 3, -1.0); 
          daxpy_(size, 2.0, d, 1, buf, 1); 

          size_t iall = 0;
          for (int j3 = i3->offset(); j3 != i3->offset()+i3->size(); ++j3) {
            for (int j2 = i2->offset(); j2 != i2->offset()+i2->size(); ++j2) {
              for (int j1 = i1->offset(); j1 != i1->offset()+i1->size(); ++j1) {
                 for (int j0 = i0->offset(); j0 != i0->offset()+i0->size(); ++j0, ++iall) {
                  en += buf[iall] * d[iall] / (eig[j0] + eig[j2] - eig[j3+nocc] - eig[j1+nocc]);
                }
              }
            }
          }

        }
      }
    }
  }
  cout << setprecision(10) << setw(20) << en << endl;

}
