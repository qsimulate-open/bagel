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
    std::unique_ptr<double[]> shalf_xx_;
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
      for (auto& i3 : virt_) {
        for (auto& i2 : active_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : active_) {
              std::vector<size_t> h = {i0.key(), i1.key(), i2.key(), i3.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              // first sort
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1(new double[r->get_size(h)]);
#if 0
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());

              // move to orthogonal basis

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
#endif
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
        // metric inverse
        const int nact = ref_->nact();
        const int nclosed = ref_->nclosed();
        const size_t dim = nact*nact;
        const size_t size = dim*dim;
        std::unique_ptr<double[]> work2(new double[size]); 
        std::unique_ptr<double[]> work(new double[std::max(size, dim*5)]);
        // TODO hardwired 0
        std::copy(ref_->rdm2(0)->data(), ref_->rdm2(0)->data()+size, work.get()); 
        sort_indices<0,2,1,3,0,1,1,1>(work, work2, nact, nact, nact, nact);
        int info;
        std::unique_ptr<double[]> eig(new double[dim]);
        dsyev_("V", "L", dim, work2.get(), dim, eig.get(), work.get(), dim*5, info);
        if (info) throw std::runtime_error("DSYEV solved");
        for (int i = 0; i != dim; ++i)
          dscal_(dim, std::pow(eig[i], -0.25), work2.get()+i*dim, 1); 
        dgemm_("N", "T", dim, dim, dim, 1.0, work2, dim, work2, dim, 0.0, work, dim); 
        shalf_xx_ = std::move(work);

        // denominator Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) * S(x2,x4; D) * S(x3, x5; D)
        // first compute Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) 
        // form f(x0,x1) <- this is not so simple...
        std::unique_ptr<double[]> work3(new double[dim*dim]);
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            std::unique_ptr<double[]> dat = this->f1_->get_block({i0.key(), i1.key()});
            for (int j1 = i1.offset(), iall = 0; j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                work3[j0-nclosed+nact*(j1-nclosed)] = dat[iall]; 
          }
        }
        // TODO hardwired 0
        dgemv_("N", size, dim, 1.0, ref_->compute_rdm3(0)->data(), size, work3.get(), 1, 0.0, work2.get(), 1);  
        // GammaF(x2,x3, x4,x5) * S(x2,x4; D) * S(x3, x5; D)
        sort_indices<0,2,1,3,0,1,1,1>(work2, work3, nact, nact, nact, nact);
        dgemm_("N", "N", dim, dim, dim, 1.0, work3, dim, shalf_xx_, dim, 0.0, work2, dim);
        dgemm_("T", "N", dim, dim, dim, 1.0, shalf_xx_, dim, work2, dim, 0.0, work3, dim);
        for (int i = 0; i != dim; ++i) eig[i] = work3[i+i*dim];
        denom_xx_ = std::move(eig);
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

