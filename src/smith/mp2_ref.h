//
// Newint - Parallel electron correlation program.
// Filename: mp2_ref.h
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


#ifndef __SRC_SMITH_MP2_REF_H 
#define __SRC_SMITH_MP2_REF_H 

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>

namespace SMITH {

// base class for Task objects
// assumes that the operation table is static (not adjustable at runtime). 
template <typename T>
class Task {
  protected:
//  std::list<std::shared_ptr<Task> > depend_;

  public:
    Task() {}; 
    ~Task() { };
    virtual void compute(std::vector<std::shared_ptr<Tensor<T> > >& o, const std::vector<IndexRange>& i) = 0;

    bool ready() const {
      bool out = true;
//    for (auto i = depend_.begin(); i != depend_.end(); ++i) out &= (*i)->ready();
      return out;
    };
};

template <typename T>
class Task0 : public Task<T> {
  public:
    Task0() : Task<T>() {};
    ~Task0() {};
    void compute(std::vector<std::shared_ptr<Tensor<T> > >& t, const std::vector<IndexRange>& i) {
      std::shared_ptr<Tensor<T> > r2_ =  t[0];
      std::shared_ptr<Tensor<T> > f2_ =  t[1];
      std::shared_ptr<Tensor<T> > t2_ =  t[2];
    };
};

template <typename T>
class MP2_Ref : public SpinFreeMethod<T> {
  protected:
//  std::list<Task> tasks_;

  public:
    MP2_Ref(std::shared_ptr<Reference> r) : SpinFreeMethod<T>(r) {
    };
    ~MP2_Ref() {};

    // get functions for private members of the base class
    IndexRange& virt() { return this->virt_; };
    IndexRange& all() { return this->all_; };
    IndexRange& closed() { return this->closed_; };
    std::shared_ptr<Reference>& ref() { return this->ref_; };;
    std::shared_ptr<Tensor<T> >& v2() { return this->v2_; };
    std::shared_ptr<Tensor<T> >& f1() { return this->f1_; };

    void mp2_iter() {

      // debug implementation of MP2 here.
      std::shared_ptr<Fock<1> > fock0(new Fock<1>(ref()->geom(), ref()->hcore()));
      std::shared_ptr<Matrix1e> den(new Matrix1e(ref()->coeff()->form_density_rhf()));
      std::shared_ptr<Fock<1> > fock1(new Fock<1>(ref()->geom(), fock0, den, ref()->schwarz()));
      Matrix1e f = *ref()->coeff() % *fock1 * *ref()->coeff();
      std::vector<double> eig;
      const int nb = ref()->nclosed() + ref()->nact() + ref()->nvirt();
      for (int i = 0; i != nb; ++i) { eig.push_back(f.element(i,i)); } // <-- of course it should be got from f1 directly :-)

      std::shared_ptr<Tensor<T> > t2 = v2()->clone();
      t2->zero();
      std::shared_ptr<Tensor<T> > r2 = t2->clone();

      //////// start iteration ///////
      for (int iter = 0; iter != 10; ++iter) {
        // setting source term
        mp2_residual(r2,t2);
        std::cout << std::setprecision(10) << std::setw(30) << mp2_energy(t2)/2 <<  "  +++" << std::endl;
        t2->daxpy(1.0, mp2_denom(r2, eig));
      }
    };

    void mp2_residual(std::shared_ptr<Tensor<T> >& r2, const std::shared_ptr<Tensor<T> > t2) {

      r2->zero();

      for (auto i3 = virt().begin(); i3 != virt().end(); ++i3) {
        for (auto i2 = closed().begin(); i2 != closed().end(); ++i2) {
          for (auto i1 = virt().begin(); i1 != virt().end(); ++i1) {
            for (auto i0 = closed().begin(); i0 != closed().end(); ++i0) {
              std::vector<size_t> h = vec(i0->key(), i1->key(), i2->key(), i3->key());
              std::vector<size_t> g = vec(i0->key(), i3->key(), i2->key(), i1->key());
              std::unique_ptr<double[]> data0 = v2()->get_block(h);
              const std::unique_ptr<double[]> data1 = v2()->get_block(g);
              sort_indices4(data1, data0, i0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 2.0, -1.0); 
              r2->put_block(h,data0);
            }
          }
        }
      }

      for (auto i3 = virt().begin(); i3 != virt().end(); ++i3) {
        for (auto i2 = closed().begin(); i2 != closed().end(); ++i2) {
          for (auto i1 = virt().begin(); i1 != virt().end(); ++i1) {
            for (auto i0 = closed().begin(); i0 != closed().end(); ++i0) {
              std::vector<size_t> ohash = vec(i0->key(), i1->key(), i2->key(), i3->key());
              std::unique_ptr<double[]> odata = r2->move_block(ohash);

              for (auto c0 = closed().begin(); c0 != closed().end(); ++c0) {
                std::vector<size_t> ihash0 = vec(c0->key(), i0->key());
                const std::unique_ptr<double[]> idata0 = f1()->get_block(ihash0);

                std::vector<size_t> ihash1 = vec(c0->key(), i1->key(), i2->key(), i3->key());
                std::vector<size_t> ihash2 = vec(c0->key(), i3->key(), i2->key(), i1->key());
                const std::unique_ptr<double[]> idata1 = t2->get_block(ihash1);
                const std::unique_ptr<double[]> idata2 = t2->get_block(ihash2);
                std::unique_ptr<double[]> idata3(new double[t2->get_size(ihash1)]); 

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


      for (auto i3 = virt().begin(); i3 != virt().end(); ++i3) {
        for (auto i2 = closed().begin(); i2 != closed().end(); ++i2) {
          for (auto i1 = virt().begin(); i1 != virt().end(); ++i1) {
            for (auto i0 = closed().begin(); i0 != closed().end(); ++i0) {
              std::vector<size_t> ohash(4); ohash[0] = i0->key(); ohash[1] = i1->key(); ohash[2] = i2->key(); ohash[3] = i3->key();
              std::unique_ptr<double[]> odata = r2->move_block(ohash);
              std::unique_ptr<double[]> idata4(new double[r2->get_size(ohash)]);

              for (auto c0 = virt().begin(); c0 != virt().end(); ++c0) {
                std::vector<size_t> ihash0 = vec(c0->key(), i0->key());
                const std::unique_ptr<double[]> idata0 = f1()->get_block(ihash0);

                std::vector<size_t> ihash1 = vec(i0->key(), c0->key(), i2->key(), i3->key());
                std::vector<size_t> ihash2 = vec(i0->key(), i3->key(), i2->key(), c0->key());
                const std::unique_ptr<double[]> idata1 = t2->get_block(ihash1);
                const std::unique_ptr<double[]> idata2 = t2->get_block(ihash2);

                std::unique_ptr<double[]> idata3(new double[t2->get_size(ihash1)]);
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
      r2 = r2->add_dagger();

    };

    std::shared_ptr<Tensor<T> > mp2_denom(std::shared_ptr<Tensor<T> > r, std::vector<double> eig) {
      std::shared_ptr<Tensor<T> > out = r->clone();
      std::vector<IndexRange> o = r->indexrange();
      assert(o.size() == 4);
      for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
        for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
          for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
            for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
              std::vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
              std::vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();
              std::unique_ptr<double[]> data0 = r->get_block(h);
              const std::unique_ptr<double[]> data1 = r->get_block(g);
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
    };


    // r is an amplitude tensor
    double mp2_energy(const std::shared_ptr<Tensor<T> > r) {
      std::shared_ptr<Tensor<T> > out = r->clone();
      std::vector<IndexRange> o = r->indexrange();
      assert(o.size() == 4);
      double en = 0.0;
      for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
        for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
          for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
            for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
              std::vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
              std::vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();
              const std::unique_ptr<double[]> data = v2()->get_block(h);
              std::unique_ptr<double[]> data0 = r->get_block(h);
              const std::unique_ptr<double[]> data1 = r->get_block(g);
              sort_indices4(data1, data0, i0->size(), i3->size(), i2->size(), i1->size(), 0, 3, 2, 1, 2.0, -1.0); 
              en += ddot_(r->get_size(h), data, 1, data0, 1);
            }
          }
        }
      }
      return en;
    };

};

}

#endif
