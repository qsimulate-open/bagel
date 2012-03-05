//
// Newint - Parallel electron correlation program.
// Filename: spinfreebase.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_SMITH_SPINFREEBASE_H
#define __SRC_SMITH_SPINFREEBASE_H

#include <src/smith/prim_op.h>
#include <src/smith/tensor.h>
#include <src/smith/moint.h>
#include <src/wfn/reference.h>

namespace SMITH {

template <typename T>
class SpinFreeMethod {
  protected:
    IndexRange virt_;
    IndexRange act_;
    IndexRange closed_;
    IndexRange all_;
    std::shared_ptr<Reference> ref_;

    std::shared_ptr<Tensor<T> > v2_;
    std::shared_ptr<Tensor<T> > f1_;

    size_t time_;

    // the diagonal denominator
    std::unique_ptr<double[]> eig_;

    // printing functions called from the solve function of a derived class
    void print_iteration() {
      std::cout << "      ---- iteration ----" << std::endl << std::endl;
      time_ = ::clock();
    };

    void print_iteration(const int i, const double en, const double err) {
      const size_t end = ::clock();
      const double tim = static_cast<double>(end - time_) / CLOCKS_PER_SEC;
      std::cout << "     " << std::setw(4) << i << std::setw(15) << std::fixed << std::setprecision(10) << en
                                                << std::setw(15) << std::fixed << std::setprecision(10) << err
                                                << std::setw(10) << std::fixed << std::setprecision(3) << tim << std::endl;
      time_ = end;
    };

    void print_iteration(const bool noconv) {
      std::cout << std::endl << "      -------------------" << std::endl;
      if (noconv) std::cout << "      *** Convergence not reached ***" << std::endl;
      std::cout << std::endl;
    };


    void update_amplitude(std::shared_ptr<Tensor<T> > t, const std::shared_ptr<Tensor<T> > r) {

      // ranks of t and r are assumed to be the same

      for (auto i3 = virt_.begin(); i3 != virt_.range().end(); ++i3) {
        for (auto i2 = closed_.begin(); i2 != closed_.range().end(); ++i2) {
          for (auto i1 = virt_.begin(); i1 != virt_.range().end(); ++i1) {
            for (auto i0 = closed_.begin(); i0 != closed_.range().end(); ++i0) {
              std::vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
              std::vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1 = r->get_block(g);

              // this is an inverse of the overlap.
              sort_indices<0,3,2,1,2,3,1,3>(data1, data0, i0->size(), i3->size(), i2->size(), i1->size()); 
              size_t iall = 0;
              for (int j3 = i3->offset(); j3 != i3->offset()+i3->size(); ++j3)
                for (int j2 = i2->offset(); j2 != i2->offset()+i2->size(); ++j2)
                  for (int j1 = i1->offset(); j1 != i1->offset()+i1->size(); ++j1)
                    for (int j0 = i0->offset(); j0 != i0->offset()+i0->size(); ++j0, ++iall)
                      data0[iall] /= (eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
              t->add_block(h,data0);
            }
          }
        }
      }
    };

  public:
    SpinFreeMethod(std::shared_ptr<Reference> r) : ref_(r) {
      const int max = 10;
      IndexRange c(r->nclosed(), max);
      IndexRange act(r->nact(), max, c.nblock(), c.size());
      IndexRange v(r->nvirt(), max, c.nblock()+act.nblock(), c.size()+act.size());
      IndexRange a(c); a.merge(act); a.merge(v);
      closed_ = c;
      act_ = act;
      virt_ = v;
      all_ = a;

      // v2 tensor.
      {
        std::vector<IndexRange> o;
        o.push_back(closed_); o.push_back(virt_);
        o.push_back(closed_); o.push_back(virt_);
        K2ext<T> v2k(ref_, o);
        v2_ = v2k.tensor();
      }
      // f1 tensor.
      {
        std::vector<IndexRange> o;
        o.push_back(all_); o.push_back(all_);
        MOFock<T> fock(ref_, o);
        f1_ = fock.tensor();
      }
    };

    IndexRange& virt() { return virt_; };
    IndexRange& all() { return all_; };
    IndexRange& closed() { return closed_; };

    std::shared_ptr<Reference>& ref() { return ref_; };;

    std::shared_ptr<Tensor<T> >& v2() { return v2_; };
    std::shared_ptr<Tensor<T> >& f1() { return f1_; };

    virtual void solve() = 0;

};

}

#endif

