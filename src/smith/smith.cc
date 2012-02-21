//
// Newint - Parallel electron correlation program.
// Filename: smith.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
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

void mp2_noniter(shared_ptr<Reference> r){
  const int max = 7;
  IndexRange closed(r->nclosed(), max);
  IndexRange virt(r->nvirt(), max, closed.nblock(), closed.size());

  vector<IndexRange> o;
  o.push_back(closed);
  o.push_back(virt);
  o.push_back(closed);
  o.push_back(virt);

  K2ext<Storage_Incore> a(r, o);
  shared_ptr<Tensor<Storage_Incore> > tensor2 = a.data();
  shared_ptr<Tensor<Storage_Incore> > tensor = tensor2->copy();

  // debug implementation of MP2 here.
  shared_ptr<Fock<1> > fock0(new Fock<1>(r->geom(), r->hcore()));
  shared_ptr<Matrix1e> den(new Matrix1e(r->coeff()->form_density_rhf()));
  shared_ptr<Fock<1> > fock1(new Fock<1>(r->geom(), fock0, den, r->schwarz()));
  Matrix1e f = *r->coeff() % *fock1 * *r->coeff();

  vector<double> eig;
  const int nb = r->nclosed() + r->nact() + r->nvirt();
  const int nocc = r->nclosed() + r->nact();
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

          sort_indices4(e, buf, i0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 0.0, -1.0); 
          daxpy_(size, 2.0, d, 1, buf, 1);

          size_t iall = 0;
          for (int j3 = i3->offset(); j3 != i3->offset()+i3->size(); ++j3) {
            for (int j2 = i2->offset(); j2 != i2->offset()+i2->size(); ++j2) {
              for (int j1 = i1->offset(); j1 != i1->offset()+i1->size(); ++j1) {
                for (int j0 = i0->offset(); j0 != i0->offset()+i0->size(); ++j0, ++iall) {
                  en += buf[iall] * d[iall] / (eig[j0] + eig[j2] - eig[j3] - eig[j1]);
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


shared_ptr<Tensor<Storage_Incore> >  mp2_denom(shared_ptr<Tensor<Storage_Incore> > r, vector<double> eig) {
  shared_ptr<Tensor<Storage_Incore> > out = r->clone();
  vector<IndexRange> o = r->indexrange();
  assert(o.size() == 4);
  for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
    for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
      for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
        for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
          vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
          vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();
          unique_ptr<double[]> data0 = r->get_block(h);
          const unique_ptr<double[]> data1 = r->get_block(g);
          sort_indices4(data1, data0, i0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 2.0/3.0, 1.0/3.0); 
          size_t iall = 0;
          for (int j3 = i3->offset(); j3 != i3->offset()+i3->size(); ++j3) {
            for (int j2 = i2->offset(); j2 != i2->offset()+i2->size(); ++j2) {
              for (int j1 = i1->offset(); j1 != i1->offset()+i1->size(); ++j1) {
                for (int j0 = i0->offset(); j0 != i0->offset()+i0->size(); ++j0, ++iall) {
                  data0[iall] /= (eig[j0] + eig[j2] - eig[j3] - eig[j1] + 0.01);
                }
              }
            }
          }
          out->put_block(h,data0);
        }
      }
    }
  }
  return out;
}


// overwrites on r (anyway this will be generated in the future)
shared_ptr<Tensor<Storage_Incore> > mp2_source(shared_ptr<Tensor<Storage_Incore> > r) {
  shared_ptr<Tensor<Storage_Incore> > out = r->clone();
  vector<IndexRange> o = r->indexrange();
  assert(o.size() == 4);
  for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
    for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
      for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
        for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
          vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
          vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();
          unique_ptr<double[]> data0 = r->get_block(h);
          const unique_ptr<double[]> data1 = r->get_block(g);
          sort_indices4(data1, data0, i0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 2.0, -1.0); 
          out->put_block(h,data0);
        }
      }
    }
  }
  return out;
}


// overwrites on r (anyway this will be generated in the future)
shared_ptr<Tensor<Storage_Incore> > mp2_add_dagger(shared_ptr<Tensor<Storage_Incore> > r) {
  shared_ptr<Tensor<Storage_Incore> > out = r->clone();
  vector<IndexRange> o = r->indexrange();
  assert(o.size() == 4);
  for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
    for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
      for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
        for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
          vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
          vector<size_t> g(4); g[0] = i2->key(); g[1] = i3->key(); g[2] = i0->key(); g[3] = i1->key();
          unique_ptr<double[]> data0 = r->get_block(h);
          const unique_ptr<double[]> data1 = r->get_block(g);
          sort_indices4(data1, data0, i2->size(), i3->size(), i0->size(), i1->size(), 2, 3, 0, 1, 1.0, 1.0); 
          out->put_block(h,data0);
        }
      }
    }
  }
  return out;
}


void mp2_iter(shared_ptr<Reference> r){ 

  const int max = 20;
  IndexRange closed(r->nclosed(), max);
  IndexRange virt(r->nvirt(), max, closed.nblock(), closed.size());
  IndexRange all(closed); all.merge(virt);

  // t2 and v2 tensors.
  vector<IndexRange> o;
  o.push_back(closed);
  o.push_back(virt);
  o.push_back(closed);
  o.push_back(virt);
  K2ext<Storage_Incore> v2k(r, o);
  const shared_ptr<Tensor<Storage_Incore> > v2 = v2k.tensor();

  // Fock
  vector<IndexRange> of;
  of.push_back(all);
  of.push_back(all);
  MOFock<Storage_Incore> a(r, of);
  shared_ptr<Tensor<Storage_Incore> > f1 = a.tensor();


  // debug implementation of MP2 here.
  shared_ptr<Fock<1> > fock0(new Fock<1>(r->geom(), r->hcore()));
  shared_ptr<Matrix1e> den(new Matrix1e(r->coeff()->form_density_rhf()));
  shared_ptr<Fock<1> > fock1(new Fock<1>(r->geom(), fock0, den, r->schwarz()));
  Matrix1e f = *r->coeff() % *fock1 * *r->coeff();
  vector<double> eig;
  const int nb = r->nclosed() + r->nact() + r->nvirt();
  for (int i = 0; i != nb; ++i) { eig.push_back(f.element(i,i)); } // <-- of course it should be got from f1 directly :-)

  shared_ptr<Tensor<Storage_Incore> > t2 = v2->clone();
  t2->zero();


  //////// start iteration ///////

  for (int iter = 0; iter != 10; ++iter) {
    // setting source term
    shared_ptr<Tensor<Storage_Incore> > r2 = v2->copy();
    r2 = mp2_source(r2);

    cout << setw(30) << r2->ddot(t2)/2 << endl;

    // test f*t2 == r2(source)

    for (auto i3 = virt.begin(); i3 != virt.end(); ++i3) {
      for (auto i2 = closed.begin(); i2 != closed.end(); ++i2) {
        for (auto i1 = virt.begin(); i1 != virt.end(); ++i1) {
          for (auto i0 = closed.begin(); i0 != closed.end(); ++i0) {
            vector<size_t> ohash(4); ohash[0] = i0->key(); ohash[1] = i1->key(); ohash[2] = i2->key(); ohash[3] = i3->key();
            unique_ptr<double[]> odata = r2->move_block(ohash);

            for (auto c0 = closed.begin(); c0 != closed.end(); ++c0) {
              vector<size_t> ihash0(2); ihash0[0] = c0->key(); ihash0[1] = i0->key();
              const unique_ptr<double[]> idata0 = f1->get_block(ihash0);

              vector<size_t> ihash1(4); ihash1[0] = c0->key(); ihash1[1] = i1->key(); ihash1[2] = i2->key(); ihash1[3] = i3->key();
              vector<size_t> ihash2(4); ihash2[0] = c0->key(); ihash2[1] = i3->key(); ihash2[2] = i2->key(); ihash2[3] = i1->key();
              const unique_ptr<double[]> idata1 = t2->get_block(ihash1);
              const unique_ptr<double[]> idata2 = t2->get_block(ihash2);
              unique_ptr<double[]> idata3(new double[t2->get_size(ihash1)]); 

              assert(t2->get_size(ihash1) == t2->get_size(ihash2));

              sort_indices4(idata2, idata3, c0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 0.0, -1.0); 
              daxpy_(t2->get_size(ihash1), 2.0, idata1, 1, idata3, 1);

              const int common = c0->size();
              const int isize0 = i0->size();
              const int isize1 = i1->size() * i2->size() * i3->size();
              dgemm_("T", "N", isize0, isize1, common, -1.0, idata0, common, idata3, common, 1.0, odata, isize0); 
            }
            r2->put_block(ohash, odata);
          }
        } 
      }
    }


    for (auto i3 = virt.begin(); i3 != virt.end(); ++i3) {
      for (auto i2 = closed.begin(); i2 != closed.end(); ++i2) {
        for (auto i1 = virt.begin(); i1 != virt.end(); ++i1) {
          for (auto i0 = closed.begin(); i0 != closed.end(); ++i0) {
            vector<size_t> ohash(4); ohash[0] = i0->key(); ohash[1] = i1->key(); ohash[2] = i2->key(); ohash[3] = i3->key();
            unique_ptr<double[]> odata = r2->move_block(ohash);
            unique_ptr<double[]> idata4(new double[r2->get_size(ohash)]);

            for (auto c0 = virt.begin(); c0 != virt.end(); ++c0) {
              vector<size_t> ihash0(2); ihash0[0] = c0->key(); ihash0[1] = i1->key();
              const unique_ptr<double[]> idata0 = f1->get_block(ihash0);

              vector<size_t> ihash1(4); ihash1[0] = i0->key(); ihash1[1] = c0->key(); ihash1[2] = i2->key(); ihash1[3] = i3->key();
              vector<size_t> ihash2(4); ihash2[0] = i0->key(); ihash2[1] = i3->key(); ihash2[2] = i2->key(); ihash2[3] = c0->key();
              const unique_ptr<double[]> idata1 = t2->get_block(ihash1);
              const unique_ptr<double[]> idata2 = t2->get_block(ihash2);

              unique_ptr<double[]> idata3(new double[t2->get_size(ihash1)]);
              sort_indices4(idata1, idata3, i0->size(), c0->size(), i2->size(), i3->size(), 1, 0, 2, 3, 0.0,  2.0);
              sort_indices4(idata2, idata3, i0->size(), i3->size(), i2->size(), c0->size(), 3, 0, 2, 1, 1.0, -1.0);

              const int common = c0->size();
              const int isize0 = i1->size();
              const int isize1 = i0->size() * i2->size() * i3->size();
              dgemm_("T", "N", isize0, isize1, common, 1.0, idata0, common, idata3, common, 0.0, idata4, isize0); 

              sort_indices4(idata4, odata, i1->size(), i0->size(), i2->size(), i3->size(), 1, 0, 2, 3, 1.0, 1.0);
            }
            r2->put_block(ohash, odata);
          }
        }
      }
    }
    r2 = mp2_add_dagger(r2);

    // this should be correct
    t2->daxpy(1.0, mp2_denom(r2, eig));

  }

}



