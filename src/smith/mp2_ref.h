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
#include <src/smith/queue.h>
#include <src/smith/mp2_ref_task.h>

namespace SMITH {

template <typename T>
class MP2_Ref : public SpinFreeMethod<T> {
  protected:
    std::shared_ptr<Queue<T> > queue_;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > r2;

  public:
    MP2_Ref(std::shared_ptr<Reference> r) : SpinFreeMethod<T>(r), queue_(new Queue<T>()) {

      t2 = this->v2_->clone();
      r2 = t2->clone();

      std::vector<std::shared_ptr<Tensor<T> > > tensor0 = vec(r2, this->f1_, t2);
      std::vector<std::shared_ptr<Tensor<T> > > tensor1 = vec(r2, this->v2_);
      std::vector<std::shared_ptr<Tensor<T> > > tensor2 = vec(r2);
      std::vector<IndexRange> index0 =  vec(this->closed_, this->virt_);

      std::shared_ptr<Task0<T> > t0(new Task0<T>(tensor0, index0));
      std::shared_ptr<Task1<T> > t1(new Task1<T>(tensor0, index0));
      std::shared_ptr<Task2<T> > t2(new Task2<T>(tensor1, index0));
      std::shared_ptr<Task3<T> > t3(new Task3<T>(tensor2, index0));
      t0->add_dep(t3);
      t1->add_dep(t3);
      t2->add_dep(t3);
      std::shared_ptr<Task4<T> > t4(new Task4<T>(tensor2, index0));
      t4->add_dep(t0);
      t4->add_dep(t1);
      t4->add_dep(t2);
      t4->add_dep(t3);

      queue_->add_task(t3);
      queue_->add_task(t4);
      queue_->add_task(t2);
      queue_->add_task(t1);
      queue_->add_task(t0);
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

      t2->zero();

      //////// start iteration ///////
      for (int iter = 0; iter != 10; ++iter) {

        queue_->initialize();
        while (!queue_->done()) queue_->next()->compute(); 

        std::cout << std::setprecision(10) << std::setw(30) << mp2_energy(t2)/2 <<  "  +++" << std::endl;
        t2->daxpy(1.0, mp2_denom(r2, eig));
      }
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
