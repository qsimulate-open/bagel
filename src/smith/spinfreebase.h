//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_SPINFREEBASE_H
#define __SRC_SMITH_SPINFREEBASE_H

#include <stddef.h>
#include <src/smith/prim_op.h>
#include <src/smith/tensor.h>
#include <src/smith/moint.h>
#include <src/util/matrix.h>
#include <src/wfn/reference.h>
#include <chrono>

namespace bagel {
namespace SMITH {

template <typename T>
class SpinFreeMethod {
  protected:
    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    IndexRange all_;
    std::shared_ptr<const Reference> ref_;

    double e0_;

    std::shared_ptr<Tensor<T> > v2_;
    std::shared_ptr<Tensor<T> > f1_;
    std::shared_ptr<Tensor<T> > rdm1_;
    std::shared_ptr<Tensor<T> > rdm2_;
    std::shared_ptr<Tensor<T> > rdm3_;
    std::shared_ptr<Tensor<T> > rdm4_;

    std::chrono::high_resolution_clock::time_point time_;

    // the diagonal denominator
    std::unique_ptr<double[]> eig_;

    // printing functions called from the solve function of a derived class
    void print_iteration() {
      std::cout << "      ---- iteration ----" << std::endl << std::endl;
      time_ = std::chrono::high_resolution_clock::now();
    };

    void print_iteration(const int i, const double en, const double err) {
      auto end = std::chrono::high_resolution_clock::now();
      const double tim = std::chrono::duration_cast<std::chrono::milliseconds>(end-time_).count() * 0.001;
      std::cout << "     " << std::setw(4) << i << std::setw(15) << std::fixed << std::setprecision(10) << en
                                                << std::setw(15) << std::fixed << std::setprecision(10) << err
                                                << std::setw(10) << std::fixed << std::setprecision(2) << tim << std::endl;
      time_ = end;
    };

    void print_iteration(const bool noconv) {
      std::cout << std::endl << "      -------------------" << std::endl;
      if (noconv) std::cout << "      *** Convergence not reached ***" << std::endl;
      std::cout << std::endl;
    };


    // E0 is defined as Trace(f(x,x), gamma(x,x))
    // For instance, E0 is 0 for MP2.
    double compute_e0() {
      if (ref_->nact() != 0 && !(static_cast<bool>(f1_) && static_cast<bool>(rdm1_)))
        throw std::logic_error("SpinFreeMethod::compute_e0 was called before f1_ or rdm1_ was computed. Strange.");
      double sum = 0.0;
      // TODO parallelize?
      for (auto& i1 : active_) {
        for (auto& i0 : active_) {
          std::vector<size_t> hash = {i0.key(), i1.key()};
          const size_t size = i0.size() * i1.size();
          std::unique_ptr<double[]> fdata = f1_->get_block(hash);
          std::unique_ptr<double[]> rdata = rdm1_->get_block(hash);
          sum += ddot_(size, fdata, 1, rdata, 1);
        }
      }
      std::cout << "    - Zeroth order energy: " << sum << std::endl;
      return sum;
    };


    // S^-1/2 for aa/xx blocks
    std::shared_ptr<Matrix> shalf_xx_;
    std::unique_ptr<double[]> denom_xx_;


    void update_amplitude(std::shared_ptr<Tensor<T> > t, const std::shared_ptr<Tensor<T> > r, const bool put = false) {

      // ranks of t and r are assumed to be the same

      // TODO should be parallelized
      for (auto& i3 : virt_) {
        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              std::vector<size_t> h = {i0.key(), i1.key(), i2.key(), i3.key()};
              std::vector<size_t> g = {i0.key(), i3.key(), i2.key(), i1.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1 = r->get_block(g);

              // this is an inverse of the overlap.
              sort_indices<0,3,2,1,2,3,1,3>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                      data0[iall] /= (eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
              if (!put) {
                t->add_block(h,data0);
              } else {
                t->put_block(h,data0);
              }
            }
          }
        }
      }
      for (auto& i2 : active_) {
        for (auto& i0 : active_) {
          // trans is the transformation matrix
          assert(shalf_xx_);
          const int nact = ref_->nact();
          std::unique_ptr<double[]> trans(new double[i0.size()*i2.size()*nact*nact]);
          for (int j2 = i2.offset(), k = 0; j2 != i2.offset()+i2.size(); ++j2)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
              std::copy_n(shalf_xx_->element_ptr(0,(j0-ref_->nclosed())+(j2-ref_->nclosed())*nact), nact*nact, trans.get()+nact*nact*k);

          for (auto& i3 : virt_) {
            for (auto& i1 : virt_) {
              std::vector<size_t> h = {i0.key(), i1.key(), i2.key(), i3.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1(new double[r->get_size(h)]);
              // sort. Active indices run faster
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // data1 is the intermediate area
              std::unique_ptr<double[]> interm(new double[i1.size()*i3.size()*nact*nact]);

              // move to orthogonal basis
              dgemm_("N", "N", nact*nact, i1.size()*i3.size(), i0.size()*i2.size(), 1.0, trans, nact*nact, data1, i0.size()*i2.size(),
                                                                                    0.0, interm, nact*nact); 

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j02 = 0; j02 != nact*nact; ++j02, ++iall)
                    interm[iall] /= e0_ - (denom_xx_[j02] + eig_[j3] + eig_[j1]);

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              dgemm_("T", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, trans, nact*nact, interm, nact*nact,
                                                                                    0.0, data1, i0.size()*i2.size()); 

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data1, data0, i0.size(), i2.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(h,data0);
              } else {
                t->put_block(h,data0);
              }
            }
          }
        }
      }
    };

  public:
    SpinFreeMethod(std::shared_ptr<const Reference> r) : ref_(r) {
      const int max = 10;
      IndexRange c(r->nclosed(), max);
      IndexRange act(r->nact(), max, c.nblock(), c.size());
      IndexRange v(r->nvirt(), max, c.nblock()+act.nblock(), c.size()+act.size());
      IndexRange a(c); a.merge(act); a.merge(v);
      closed_ = c;
      active_ = act;
      virt_ = v;
      all_ = a;

      // v2 tensor.
      {
        IndexRange occ(closed_);
        occ.merge(active_);
        IndexRange virt(active_);
        virt.merge(virt_);

        std::vector<IndexRange> o = {occ, virt, occ, virt};
        K2ext<T> v2k(ref_, o);
        v2_ = v2k.tensor();
      }
      // f1 tensor.
      {
        std::vector<IndexRange> o = {all_, all_};
        MOFock<T> fock(ref_, o);
        f1_ = fock.tensor();
      }
      // rdms
      if (!ref_->rdm1().empty()) {
        std::vector<IndexRange> o = {active_, active_};
        rdm1_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false)); 
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            std::vector<size_t> hash = {i0.key(), i1.key()};
            const size_t size = i0.size() * i1.size();
            std::unique_ptr<double[]> data(new double[size]); 
            int iall = 0;
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                // TODO for the time being we hardwire "0" here (but this should be fixed)
                data[iall] = ref_->rdm1(0)->element({j0, j1});
            rdm1_->put_block(hash, data);
          }
        }
      }
      if (!ref_->rdm2().empty()) {
        std::vector<IndexRange> o = {active_, active_, active_, active_};
        rdm2_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false)); 
        for (auto& i3 : active_) {
          for (auto& i2 : active_) {
            for (auto& i1 : active_) {
              for (auto& i0 : active_) {
                std::vector<size_t> hash = {i0.key(), i1.key(), i2.key(), i3.key()};
                const size_t size = i0.size() * i1.size() * i2.size() * i3.size();
                std::unique_ptr<double[]> data(new double[size]); 
                int iall = 0;
                for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                        // TODO for the time being we hardwire "0" here (but this should be fixed)
                        data[iall] = ref_->rdm2(0)->element({j0, j1, j2, j3});
                 rdm2_->put_block(hash, data);
              }
            }
          }
        }
      }
      // generic function??
      if (!ref_->rdm1().empty() && !ref_->rdm2().empty()) {
        std::vector<IndexRange> o = {active_, active_, active_, active_, active_, active_};
        rdm3_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false));
        // TODO for the time being we hardwire "0" here (but this should be fixed)
        std::shared_ptr<RDM<3> > rdm3source = ref_->compute_rdm3(0);
        for (auto& i5 : active_) {
          for (auto& i4 : active_) {
            for (auto& i3 : active_) {
              for (auto& i2 : active_) {
                for (auto& i1 : active_) {
                  for (auto& i0 : active_) {
                    std::vector<size_t> hash = {i0.key(), i1.key(), i2.key(), i3.key(), i4.key(), i5.key()};
                    const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size();
                    std::unique_ptr<double[]> data(new double[size]); 
                    int iall = 0;
                    for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                      for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                        for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                          for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                                data[iall] = rdm3source->element({j0, j1, j2, j3, j4, j5});
                    rdm3_->put_block(hash, data);
                  }
                }
              }
            }
          }
        }
      }

      // aa/xx blocks 
      if (!ref_->rdm2().empty()) {
        // metric half inverse (S^-1/2)
        const int nact = ref_->nact();
        const int nclosed = ref_->nclosed();
        const size_t dim = nact*nact;
        const size_t size = dim*dim;
        std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
        {
          // TODO hardwired 0
          shalf_xx_ = std::shared_ptr<Matrix>(new Matrix(dim, dim));
          std::copy_n(ref_->rdm2(0)->data(), size, shalf_xx_->data()); 
          sort_indices<0,2,1,3,0,1,1,1>(shalf_xx_->data(), work2->data(), nact, nact, nact, nact);
          work2->inverse_half(1.0e-9);
          *shalf_xx_ = *work2;
//#define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
          Matrix r(dim, dim);
          Matrix s(dim, dim);
          std::copy_n(ref_->rdm2(0)->data(), size, r.data());
          sort_indices<0,2,1,3,0,1,1,1>(r.data(), s.data(), nact, nact, nact, nact);
          (s * *shalf_xx_ * *shalf_xx_).print("unit");
#endif
        }

        // denominator Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) * T(x2,x4; D) * T(x3, x5; D)
        // first compute Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) // TODO this should be computed directly maybe 
        // form f(x0,x1) <- this is not so simple...
        std::shared_ptr<Matrix> fockact(new Matrix(nact, nact));
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            std::unique_ptr<double[]> dat = this->f1_->get_block({i0.key(), i1.key()});
            for (int j1 = i1.offset(), iall = 0; j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                fockact->element(j0-nclosed, j1-nclosed) = dat[iall];
          }
        }
        // TODO hardwired 0
        std::shared_ptr<RDM<3> > rdm3 = ref_->compute_rdm3(0);
        dgemv_("N", size, dim, 1.0, rdm3->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);
        // GammaF(x2,x3, x4,x5) * T(x2,x4; D) * T(x3, x5; D)
        std::shared_ptr<Matrix> work4(new Matrix(dim, dim));
        sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work4->data(), nact, nact, nact, nact);

        denom_xx_ = (*shalf_xx_ % *work4 * *shalf_xx_).diag();

#ifdef LOCAL_DEBUG
for (int i = 0; i != dim; ++i) std::cout << denom_xx_[i] << std::endl;
#endif
      }

      // set e0
      e0_ = compute_e0();
    };

    IndexRange& virt() { return virt_; };
    IndexRange& all() { return all_; };
    IndexRange& closed() { return closed_; };

    std::shared_ptr<const Reference>& ref() { return ref_; };;

    double e0() const { return e0_; };

    virtual void solve() = 0;

};

}
}

#endif

