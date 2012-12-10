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
    }

    void print_iteration(const int i, const double en, const double err) {
      auto end = std::chrono::high_resolution_clock::now();
      const double tim = std::chrono::duration_cast<std::chrono::milliseconds>(end-time_).count() * 0.001;
      std::cout << "     " << std::setw(4) << i << std::setw(15) << std::fixed << std::setprecision(10) << en
                                                << std::setw(15) << std::fixed << std::setprecision(10) << err
                                                << std::setw(10) << std::fixed << std::setprecision(2) << tim << std::endl;
      time_ = end;
    }

    void print_iteration(const bool noconv) {
      std::cout << std::endl << "      -------------------" << std::endl;
      if (noconv) std::cout << "      *** Convergence not reached ***" << std::endl;
      std::cout << std::endl;
    }


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
      std::cout << "    - Zeroth order energy: " << std::setw(20) << std::setprecision(10) << sum << std::endl;
      return sum;
    }


    // S^-1/2 for aa/xx blocks (overlap is 2rdm)
    std::shared_ptr<Matrix> shalf_xx_;
    std::unique_ptr<double[]> denom_xx_;

    // S^-1/2 for aa/cx blocks (overlap is 1rdm)
    std::shared_ptr<Matrix> shalf_x_;
    std::unique_ptr<double[]> denom_x_;

    // S^-1/2 for ax/cc blocks
    std::shared_ptr<Matrix> shalf_h_;
    std::unique_ptr<double[]> denom_h_;

    // S^-1/2 for xx/cc blocks
    std::shared_ptr<Matrix> shalf_hh_;
    std::unique_ptr<double[]> denom_hh_;

    // S^-1/2 for ax/cx and xa/cx blocks (4*2rdm size)
    std::shared_ptr<Matrix> shalf_xh_;
    std::unique_ptr<double[]> denom_xh_;

    // S^-1/2 for xa/xx blocks (3rdm size)
    std::shared_ptr<Matrix> shalf_xhh_;
    std::unique_ptr<double[]> denom_xhh_;

    // S^-1/2 for xx/cx blocks (3rdm size)
    std::shared_ptr<Matrix> shalf_xxh_;
    std::unique_ptr<double[]> denom_xxh_;

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
              const std::unique_ptr<double[]> data1 = r->get_block(g);

              // this is an inverse of the overlap.
              // prefactor of 0.25 included here
              sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                      // note that e0 is cancelled by another term
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
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i0.size()*i2.size()*nact*nact]);
          for (int j2 = i2.offset(), k = 0; j2 != i2.offset()+i2.size(); ++j2)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
              std::copy_n(shalf_xx_->element_ptr(0,(j0-nclo)+(j2-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

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
              // intermediate area
              std::unique_ptr<double[]> interm(new double[i1.size()*i3.size()*nact*nact]);

              // move to orthogonal basis
              dgemm_("N", "N", nact*nact, i1.size()*i3.size(), i0.size()*i2.size(), 1.0, transp, nact*nact, data1, i0.size()*i2.size(),
                                                                                    0.0, interm, nact*nact);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j02 = 0; j02 != nact*nact; ++j02, ++iall)
                    interm[iall] /= e0_ - (denom_xx_[j02] + eig_[j3] + eig_[j1]);

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              dgemm_("T", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, transp, nact*nact, interm, nact*nact,
                                                                                    0.0, data0,  i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(h,data1);
              } else {
                t->put_block(h,data1);
              }
            }
          }
        }
      }
      for (auto& i0 : active_) {
        // trans is the transformation matrix
        assert(shalf_x_);
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        std::unique_ptr<double[]> transp(new double[i0.size()*nact]);
        for (int j0 = i0.offset(), k = 0; j0 != i0.offset()+i0.size(); ++j0, ++k)
          std::copy_n(shalf_x_->element_ptr(0,j0-nclo), nact, transp.get()+nact*k);

        for (auto& i3 : virt_) {
          for (auto& i2 : closed_) {
            for (auto& i1 : virt_) {
              std::vector<size_t> h = {i2.key(), i3.key(), i0.key(), i1.key()};
              std::vector<size_t> g = {i2.key(), i1.key(), i0.key(), i3.key()};
              if (!r->get_size(h)) continue;
              assert(r->get_size(g));
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data2(new double[r->get_size(h)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
              const std::unique_ptr<double[]> data1 = r->get_block(g);
              sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());
              std::unique_ptr<double[]> interm(new double[i1.size()*i2.size()*i3.size()*nact]);

              // move to orthogonal basis
              dgemm_("N", "N", nact, i1.size()*i2.size()*i3.size(), i0.size(), 1.0, transp, nact, data2, i0.size(),
                                                                               0.0, interm, nact);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = 0; j0 != nact; ++j0, ++iall)
                      interm[iall] /= e0_ - (denom_x_[j0] + eig_[j3] - eig_[j2] + eig_[j1]);

              // move back to non-orthogonal basis
              dgemm_("T", "N", i0.size(), i1.size()*i2.size()*i3.size(), nact, 1.0, transp, nact, interm, nact,
                                                                               0.0, data2,  i0.size());

              std::vector<size_t> hout = {i0.key(), i1.key(), i2.key(), i3.key()};
              if (!put) {
                t->add_block(hout,data2);
              } else {
                t->put_block(hout,data2);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        // trans is the transformation matrix
        assert(shalf_h_);
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        std::unique_ptr<double[]> transp(new double[i3.size()*nact]);
        for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3, ++k)
          std::copy_n(shalf_h_->element_ptr(0,j3-nclo), nact, transp.get()+nact*k);

        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              std::vector<size_t> h = {i2.key(), i3.key(), i0.key(), i1.key()};
              std::vector<size_t> g = {i0.key(), i3.key(), i2.key(), i1.key()};
              if (!r->get_size(h)) continue;
              assert(r->get_size(g));
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data2(new double[r->get_size(h)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
              const std::unique_ptr<double[]> data1 = r->get_block(g);
              sort_indices<0,3,2,1,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
              std::unique_ptr<double[]> interm(new double[i0.size()*i1.size()*i2.size()*nact]);

              // move to orthogonal basis
              dgemm_("N", "T", i0.size()*i1.size()*i2.size(), nact, i3.size(), 1.0, data2, i0.size()*i1.size()*i2.size(), transp, nact,
                                                                               0.0, interm, i0.size()*i1.size()*i2.size());

              size_t iall = 0;
              for (int j3 = 0; j3 != nact; ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                      interm[iall] /= e0_ - (denom_h_[j3] - eig_[j2] + eig_[j1] - eig_[j0]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i0.size()*i1.size()*i2.size(), i3.size(), nact, 1.0, interm, i0.size()*i1.size()*i2.size(), transp, nact,
                                                                               0.0, data2,  i0.size()*i1.size()*i2.size());

              std::vector<size_t> hout = {i0.key(), i1.key(), i2.key(), i3.key()};
              if (!put) {
                t->add_block(hout,data2);
              } else {
                t->put_block(hout,data2);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          assert(shalf_hh_);
          const int nact = ref_->nact();
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i1.size()*i3.size()*nact*nact]);
          for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++k)
              std::copy_n(shalf_hh_->element_ptr(0,(j1-nclo)+(j3-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

          for (auto& i2 : closed_) {
            for (auto& i0 : closed_) {
              std::vector<size_t> h = {i0.key(), i1.key(), i2.key(), i3.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1(new double[r->get_size(h)]);
              // sort. Active indices run slower
              sort_indices<0,3,2,1,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
              // intermediate area
              std::unique_ptr<double[]> interm(new double[i0.size()*i2.size()*nact*nact]);

              // move to orthogonal basis
              dgemm_("N", "T", i0.size()*i2.size(), nact*nact, i1.size()*i3.size(), 1.0, data1, i0.size()*i2.size(), transp, nact*nact,
                                                                                    0.0, interm, i0.size()*i2.size());

              size_t iall = 0;
              for (int j13 = 0; j13 != nact*nact; ++j13)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= e0_ - (denom_hh_[j13] - eig_[j2] - eig_[j0]);

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              dgemm_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, interm, i0.size()*i2.size(), transp, nact*nact,
                                                                                    0.0, data0,  i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(h,data1);
              } else {
                t->put_block(h,data1);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          assert(shalf_xh_);
          const int nact = ref_->nact();
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i2.size()*i3.size()*nact*nact*4]);
          for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++k) {
              std::copy_n(shalf_xh_->element_ptr(0,             (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*k);
              std::copy_n(shalf_xh_->element_ptr(0, nact*nact + (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*(k+i2.size()*i3.size()));
            }

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {

              std::vector<size_t> hin = {i2.key(), i3.key(), i0.key(), i1.key()};
              std::vector<size_t> gin = {i0.key(), i3.key(), i2.key(), i1.key()};
              // if this block is not included in the current wave function, skip it
              if (!r->get_size(hin)) continue;
              assert(r->get_size(gin) == r->get_size(hin));
              const size_t blocksize = r->get_size(hin);
              std::unique_ptr<double[]> data0 = r->get_block(hin);
              std::unique_ptr<double[]> data1 = r->get_block(gin);

              std::unique_ptr<double[]> data2(new double[blocksize*2]);
              // sort. Active indices run slower
              sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
              sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
              // intermediate area
              std::unique_ptr<double[]> interm(new double[i0.size()*i1.size()*nact*nact*2]);

              // move to orthogonal basis
              dgemm_("N", "T", i0.size()*i1.size(), nact*nact*2, i2.size()*i3.size()*2, 1.0, data2,  i0.size()*i1.size(), transp, nact*nact*2,
                                                                                        0.0, interm, i0.size()*i1.size());

              size_t iall = 0;
              for (int j23 = 0; j23 != nact*nact*2; ++j23)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    interm[iall] /= e0_ - (denom_xh_[j23] + eig_[j1] - eig_[j0]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i0.size()*i1.size(), i2.size()*i3.size()*2, nact*nact*2, 1.0, interm, i0.size()*i1.size(), transp, nact*nact*2,
                                                                                        0.0, data2,  i0.size()*i1.size());

              // sort back to the original order
              std::copy_n(data2.get(), blocksize, data0.get());
              sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksize, data1.get(), i0.size(), i1.size(), i2.size(), i3.size());
              std::vector<size_t> h = {i0.key(), i1.key(), i2.key(), i3.key()};
              std::vector<size_t> g = {i2.key(), i1.key(), i0.key(), i3.key()};
              if (!put) {
                t->add_block(h,data0);
                t->add_block(g,data1);
              } else {
                t->put_block(h,data0);
                t->put_block(g,data1);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          for (auto& i0 : active_) {
            assert(shalf_xhh_);
            const int nact = ref_->nact();
            const int nclo = ref_->nclosed();
            std::unique_ptr<double[]> transp(new double[i0.size()*i2.size()*i3.size()*nact*nact*nact]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                  std::copy_n(shalf_xhh_->element_ptr(0,j0-nclo+nact*(j2-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

            for (auto& i1 : virt_) {
              std::vector<size_t> h = {i2.key(), i3.key(), i0.key(), i1.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1(new double[r->get_size(h)]);
              // sort. Active indices run slower
              sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
              // intermediate area
              std::unique_ptr<double[]> interm(new double[i1.size()*nact*nact*nact]);

              // move to orthogonal basis
              dgemm_("N", "T", i1.size(), nact*nact*nact, i0.size()*i2.size()*i3.size(), 1.0, data1,  i1.size(), transp, nact*nact*nact,
                                                                                         0.0, interm, i1.size());

              size_t iall = 0;
              for (int j123 = 0; j123 != nact*nact*nact; ++j123)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++iall)
                  interm[iall] /= e0_ - (denom_xhh_[j123] + eig_[j1]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i1.size(), i0.size()*i2.size()*i3.size(), nact*nact*nact, 1.0, interm, i1.size(), transp, nact*nact*nact,
                                                                                         0.0, data0,  i1.size());

              // sort back to the original order
              sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i1.size(), i0.size(), i2.size(), i3.size());
              std::vector<size_t> hout = {i0.key(), i1.key(), i2.key(), i3.key()};
              if (!put) {
                t->add_block(hout,data1);
              } else {
                t->put_block(hout,data1);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            assert(shalf_xhh_);
            const int nact = ref_->nact();
            const int nclo = ref_->nclosed();
            std::unique_ptr<double[]> transp(new double[i0.size()*i1.size()*i3.size()*nact*nact*nact]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                  std::copy_n(shalf_xxh_->element_ptr(0,j0-nclo+nact*(j1-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

            for (auto& i2 : closed_) {
              std::vector<size_t> h = {i2.key(), i3.key(), i0.key(), i1.key()};

              // if this block is not included in the current wave function, skip it
              if (!r->get_size(h)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(h);
              std::unique_ptr<double[]> data1(new double[r->get_size(h)]);
              // sort. Active indices run slower
              sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
              // intermediate area
              std::unique_ptr<double[]> interm(new double[i2.size()*nact*nact*nact]);

              // move to orthogonal basis
              dgemm_("N", "T", i2.size(), nact*nact*nact, i0.size()*i1.size()*i3.size(), 1.0, data1,  i2.size(), transp, nact*nact*nact,
                                                                                         0.0, interm, i2.size());

              size_t iall = 0;
              for (int j013 = 0; j013 != nact*nact*nact; ++j013)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++iall)
//                interm[iall] /= e0_ - (- eig_[j2]);
                  interm[iall] /= e0_ - (denom_xxh_[j013] - eig_[j2]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i2.size(), i0.size()*i1.size()*i3.size(), nact*nact*nact, 1.0, interm, i2.size(), transp, nact*nact*nact,
                                                                                         0.0, data0,  i2.size());

              // sort back to the original order
              sort_indices<1,2,0,3,0,1,1,1>(data0, data1, i2.size(), i0.size(), i1.size(), i3.size());
              std::vector<size_t> hout = {i0.key(), i1.key(), i2.key(), i3.key()};
              if (!put) {
                t->add_block(hout,data1);
              } else {
                t->put_block(hout,data1);
              }
            }
          }
        }
      }
    }

  public:
    SpinFreeMethod(std::shared_ptr<const Reference> r) : ref_(r) {
      const int max = 10;
      IndexRange c(r->nclosed(), max);
      IndexRange act(r->nact(), max, c.nblock(), c.size());
      IndexRange v(r->nvirt(), max, c.nblock()+act.nblock(), c.size()+act.size());
      IndexRange a(c); a.merge(act); a.merge(v);
      assert(c.size() == r->nclosed());
      assert(act.size() == r->nact());
      assert(v.size() == r->nvirt());
      closed_ = c;
      active_ = act;
      virt_ = v;
      all_ = a;

      // f1 tensor.
      std::shared_ptr<const Coeff> coeff;
      {
        std::vector<IndexRange> o = {all_, all_};
        MOFock<T> fock(ref_, o);
        f1_ = fock.tensor();
        // canonical orbitals within closed and virtual subspaces
        coeff = fock.coeff();
      }

      // v2 tensor.
      {
        IndexRange occ(closed_);
        occ.merge(active_);
        IndexRange virt(active_);
        virt.merge(virt_);

        std::vector<IndexRange> o = {occ, virt, occ, virt};
        K2ext<T> v2k(ref_, coeff, o);
        v2_ = v2k.tensor();
      }
      // rdms
      if (!ref_->rdm1().empty()) {
        std::vector<IndexRange> o = {active_, active_};
        rdm1_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false));
        const int nclo = ref_->nclosed();
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            std::vector<size_t> hash = {i0.key(), i1.key()};
            const size_t size = i0.size() * i1.size();
            std::unique_ptr<double[]> data(new double[size]);
            int iall = 0;
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                // TODO for the time being we hardwire "0" here (but this should be fixed)
                data[iall] = ref_->rdm1(0)->element(j0-nclo, j1-nclo);
            rdm1_->put_block(hash, data);
          }
        }
      }
      if (!ref_->rdm2().empty()) {
        std::vector<IndexRange> o = {active_, active_, active_, active_};
        rdm2_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false));
        const int nclo = ref_->nclosed();
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
                        data[iall] = ref_->rdm2(0)->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo);
                 rdm2_->put_block(hash, data);
              }
            }
          }
        }
      }
      // TODO generic function??
      if (!ref_->rdm1().empty() && !ref_->rdm2().empty()) {
        {
          std::vector<IndexRange> o = {active_, active_, active_, active_, active_, active_};
          rdm3_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(o, false));
          std::vector<IndexRange> p = {active_, active_, active_, active_, active_, active_, active_, active_};
          rdm4_ = std::shared_ptr<Tensor<T> >(new Tensor<T>(p, false));
        }

        // TODO for the time being we hardwire "0" here (but this should be fixed)
        std::shared_ptr<RDM<3> > rdm3source;
        std::shared_ptr<RDM<4> > rdm4source;
        std::tie(rdm3source, rdm4source) = ref_->compute_rdm34(0);

        const int nclo = ref_->nclosed();
        for (auto& i5 : active_)
          for (auto& i4 : active_)
            for (auto& i3 : active_)
              for (auto& i2 : active_)
                for (auto& i1 : active_)
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
                                data[iall] = rdm3source->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo);
                    rdm3_->put_block(hash, data);
                  }
        // TODO there should be a better way of doing this!!!
        for (auto& i7 : active_)
          for (auto& i6 : active_)
            for (auto& i5 : active_)
              for (auto& i4 : active_)
                for (auto& i3 : active_)
                  for (auto& i2 : active_)
                    for (auto& i1 : active_)
                      for (auto& i0 : active_) {
                        std::vector<size_t> hash = {i0.key(), i1.key(), i2.key(), i3.key(), i4.key(), i5.key(), i6.key(), i7.key()};
                        const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size() * i6.size() * i7.size();
                        std::unique_ptr<double[]> data(new double[size]);
                        int iall = 0;
                        for (int j7 = i7.offset(); j7 != i7.offset()+i7.size(); ++j7)
                          for (int j6 = i6.offset(); j6 != i6.offset()+i6.size(); ++j6)
                            for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                              for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                                for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                                        data[iall] = rdm4source->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo, j6-nclo, j7-nclo);
                    rdm4_->put_block(hash, data);
                  }

        const int nact = ref_->nact();
        std::shared_ptr<Matrix> fockact(new Matrix(nact, nact));
        for (auto& i1 : active_)
          for (auto& i0 : active_)
            fockact->copy_block(i0.offset()-nclo, i1.offset()-nclo, i0.size(), i1.size(), this->f1_->get_block({i0.key(), i1.key()}));

        // TODO to be cleaned up
        std::shared_ptr<Matrix> rdm1mat(new Matrix(nact, nact));
        std::shared_ptr<Matrix> rdm2mat(new Matrix(nact*nact, nact*nact));
        // TODO hardwired 0
        std::copy_n(ref_->rdm1(0)->data(), rdm1mat->size(), rdm1mat->data());
        std::copy_n(ref_->rdm2(0)->data(), rdm2mat->size(), rdm2mat->data());
        std::shared_ptr<RDM<2> > rdm2(new RDM<2>(*ref_->rdm2(0)));

        // aa/xx blocks
        // metric half inverse (S^-1/2)
        {
          const size_t dim = nact*nact;
          const size_t size = dim*dim;
          shalf_xx_ = std::shared_ptr<Matrix>(new Matrix(dim, dim));
          Matrix tmp = *rdm2mat;
          sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), shalf_xx_->data(), nact, nact, nact, nact);
          shalf_xx_->inverse_half(1.0e-9);

          // denominator Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) * T(x2,x4; D) * T(x3, x5; D)
          // first compute Gamma(x0,x1, x2,x3, x4,x5) * f(x0,x1) // TODO this should be computed directly maybe
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, rdm3source->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);

          // GammaF(x2,x3, x4,x5) * T(x2,x4; D) * T(x3, x5; D)
          std::shared_ptr<Matrix> work4 = work2->clone();
          sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work4->data(), nact, nact, nact, nact);
          Matrix fss = *shalf_xx_ % *work4 * *shalf_xx_;
          denom_xx_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_xx_.get());
          *shalf_xx_ = fss % *shalf_xx_;
        }

        // aa/cx blocks
        {
          const size_t dim = nact;
          const size_t size = dim*dim;
          shalf_x_ = std::shared_ptr<Matrix>(new Matrix(*rdm1mat));
          shalf_x_->inverse_half(1.0e-9);

          // denominator Gamma(x0,x1, x2,x3) * f(x0,x1) * T(x2,D) * T(x3,D)
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, rdm2->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);

          Matrix fss = *shalf_x_ % *work2 * *shalf_x_;
          denom_x_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_x_.get());
          *shalf_x_ = fss % *shalf_x_;
        }

        // ax/cc blocks
        {
          const size_t dim = nact;
          const size_t size = dim*dim;
          shalf_h_ = std::shared_ptr<Matrix>(new Matrix(*rdm1mat));
          shalf_h_->scale(-1.0);
          shalf_h_->add_diag(2.0); //.. making hole 1RDM
          std::shared_ptr<const Matrix> h1(new Matrix(*shalf_h_));
          shalf_h_->inverse_half(1.0e-9);

          // denominator hole(x0,x1, x2,x3) * f(x0,x1) * T(x2,D) * T(x3,D)
          std::shared_ptr<RDM<2> > hole(new RDM<2>(*rdm2));
          hole->scale(-1.0);
          for (int i = 0; i != nact; ++i) {
            for (int j = 0; j != nact; ++j) {
              for (int k = 0; k != nact; ++k) {
                // see Celani eq. A7
                hole->element(k, i, i, j) += h1->element(k,j);
                hole->element(i, k, j, i) -= rdm1mat->element(k,j);
                hole->element(i, i, k, j) += 2.0*rdm1mat->element(k,j);
              }
            }
          }
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, hole->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);

          Matrix fss = *shalf_h_ % *work2 * *shalf_h_;
          denom_h_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_h_.get());
          *shalf_h_ = fss % *shalf_h_;
        }
        // xx/cc blocks
        {
          // overlap (2rdm size)
          const size_t dim = nact*nact;
          const size_t size = dim*dim;
          std::shared_ptr<RDM<2> > ovl(new RDM<2>(*rdm2));
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i1 = 0; i1 != nact; ++i1)
              for (int i3 = 0; i3 != nact; ++i3)
                for (int i0 = 0; i0 != nact; ++i0) {
                  if (i1 == i3) ovl->element(i0, i3, i1, i2) +=      rdm1mat->element(i0, i2);
                  if (i1 == i2) ovl->element(i0, i3, i1, i2) += -2.0*rdm1mat->element(i0, i3);
                  if (i0 == i3) ovl->element(i0, i3, i1, i2) += -2.0*rdm1mat->element(i1, i2);
                  if (i0 == i2) ovl->element(i0, i3, i1, i2) +=      rdm1mat->element(i1, i3);
                  if (i0 == i3 && i1 == i2) ovl->element(i0, i3, i1, i2) +=  4.0;
                  if (i0 == i2 && i1 == i3) ovl->element(i0, i3, i1, i2) += -2.0;
                }
          shalf_hh_ = std::shared_ptr<Matrix>(new Matrix(dim, dim));
          sort_indices<0,2,1,3,0,1,1,1>(ovl->data(), shalf_hh_->data(), nact, nact, nact, nact);
          shalf_hh_->inverse_half(1.0e-9);
          // denominator
          std::shared_ptr<RDM<3> > rdm3(new RDM<3>(*rdm3source));
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i3 = 0; i3 != nact; ++i3)
              for (int i4 = 0; i4 != nact; ++i4)
                for (int i1 = 0; i1 != nact; ++i1)
                  for (int i5 = 0; i5 != nact; ++i5)
                    for (int i0 = 0; i0 != nact; ++i0) {
                      if (i3 == i5)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i1, i4, i0, i2);
                      if (i3 == i4)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i0, i5, i1, i2);
                      if (i1 == i5)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i0, i4, i3, i2);
                      if (i1 == i5 && i3 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i0, i2);
                      if (i1 == i4)                         rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm2->element(i0, i5, i3, i2);
                      if (i1 == i4 && i3 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i0, i2);
                      if (i1 == i2)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i0, i5, i3, i4);
                      if (i1 == i2 && i3 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i0, i4);
                      if (i1 == i2 && i3 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i0, i5);
                      if (i0 == i5)                         rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm2->element(i1, i4, i3, i2);
                      if (i0 == i5 && i3 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i1, i2);
                      if (i1 == i2 && i0 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i3, i4);
                      if (i1 == i2 && i0 == i5 && i3 == i4) rdm3->element(i0, i5, i1, i4, i3, i2) += 4.0;
                      if (i0 == i4)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i1, i5, i3, i2);
                      if (i0 == i4 && i3 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i1, i2);
                      if (i1 == i2 && i0 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i3, i5);
                      if (i1 == i2 && i0 == i4 && i3 == i5) rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0;
                      if (i0 == i2)                         rdm3->element(i0, i5, i1, i4, i3, i2) += rdm2->element(i3, i5, i1, i4);
                      if (i0 == i2 && i3 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i1, i4);
                      if (i0 == i2 && i3 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i1, i5);
                      if (i1 == i5 && i0 == i2)             rdm3->element(i0, i5, i1, i4, i3, i2) += rdm1mat->element(i3, i4);
                      if (i0 == i2 && i1 == i5 && i3 == i4) rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0;
                      if (i0 == i2 && i1 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i3, i5);
                      if (i1 == i4 && i3 == i5 && i0 == i2) rdm3->element(i0, i5, i1, i4, i3, i2) += 4.0;
                      if (i0 == i5 && i1 == i4)             rdm3->element(i0, i5, i1, i4, i3, i2) += 4.0* rdm1mat->element(i3, i2);
                      if (i0 == i4 && i1 == i5)             rdm3->element(i0, i5, i1, i4, i3, i2) += -2.0* rdm1mat->element(i3, i2);
                    }
          // denominator rdm3(x0,x1, x2,x3, x4,x5) * f(x0,x1) * T(x2,x4; D) * T(x3, x5; D)
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, rdm3->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);

          // GammaF(x2,x3, x4,x5) * T(x2,x4; D) * T(x3, x5; D)
          std::shared_ptr<Matrix> work4 = work2->clone();
          sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work4->data(), nact, nact, nact, nact);
          Matrix fss = *shalf_hh_ % *work4 * *shalf_hh_;
          denom_hh_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_hh_.get());
          *shalf_hh_ = fss % *shalf_hh_;
        }

        // ax/xi & xa/xi cases (they are treated at once)
        {
          // overlap (2rdm size)
          const size_t dim = nact*nact;
          const size_t size = dim*dim;
          shalf_xh_ = std::shared_ptr<Matrix>(new Matrix(dim*2, dim*2));
          std::shared_ptr<RDM<2> > ovl1(new RDM<2>(*rdm2));
          ovl1->scale(-1.0);
          std::shared_ptr<RDM<2> > ovl4(new RDM<2>(*rdm2));
          for (int i = 0; i != nact; ++i)
            for (int j = 0; j != nact; ++j)
              for (int k = 0; k != nact; ++k) {
                ovl1->element(i, i, k, j) += 2.0 * rdm1mat->element(k, j);
                ovl4->element(k, i, i, j) += rdm1mat->element(k, j);
              }
          std::shared_ptr<Matrix> work(new Matrix(dim, dim));
          sort_indices<3,0,2,1,0,1,1,1>(ovl1->data(), work->data(), nact, nact, nact, nact);
          shalf_xh_->add_block(dim, dim, dim, dim, work);
          sort_indices<0,1,3,2,0,1,2,1>(ovl4->data(), work->data(), nact, nact, nact, nact);
          shalf_xh_->add_block(0, 0, dim, dim, work);
          work->scale(-0.5);
          shalf_xh_->add_block(dim, 0, dim, dim, work);
          shalf_xh_->add_block(0, dim, dim, dim, work);
          shalf_xh_->inverse_half(1.0e-9);
          // denominator
          std::shared_ptr<RDM<3> > den0(new RDM<3>(*rdm3source));
          den0->scale(-1.0);
          std::shared_ptr<RDM<3> > den3(new RDM<3>(*rdm3source));
          for (int i5 = 0; i5 != nact; ++i5)
            for (int i4 = 0; i4 != nact; ++i4)
              for (int i3 = 0; i3 != nact; ++i3)
                for (int i2 = 0; i2 != nact; ++i2)
                  for (int i1 = 0; i1 != nact; ++i1)
                    for (int i0 = 0; i0 != nact; ++i0) {
                      if (i3 == i4)             den3->element(i0, i1, i4, i5, i2, i3) += rdm2->element(i0, i1, i2, i5);
                      if (i1 == i2)             den3->element(i0, i1, i4, i5, i2, i3) += rdm2->element(i0, i3, i4, i5);
                      if (i1 == i2 && i3 == i4) den3->element(i0, i1, i4, i5, i2, i3) += rdm1mat->element(i0, i5);
                      if (i1 == i4)             den3->element(i0, i1, i4, i5, i2, i3) += rdm2->element(i2, i3, i0, i5);
                      if (i3 == i4)             den0->element(i4, i1, i0, i5, i2, i3) += -1.0 * rdm2->element(i2, i1, i0, i5);
                      if (i1 == i2)             den0->element(i4, i1, i0, i5, i2, i3) += -1.0 * rdm2->element(i4, i3, i0, i5);
                      if (i3 == i4 && i1 == i2) den0->element(i4, i1, i0, i5, i2, i3) +=  2.0 * rdm1mat->element(i0, i5);
                      if (i1 == i4)             den0->element(i4, i1, i0, i5, i2, i3) +=  2.0 * rdm2->element(i2, i3, i0, i5);
                    }
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, den0->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);
          sort_indices<3,0,2,1,0,1,1,1>(work2->data(), work->data(), nact, nact, nact, nact);
          Matrix num(dim*2, dim*2);
          num.add_block(dim, dim, dim, dim, work);
          dgemv_("N", size, nact*nact, 1.0, den3->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);
          sort_indices<0,1,3,2,0,1,2,1>(work2->data(), work->data(), nact, nact, nact, nact);
          num.add_block(0, 0, dim, dim, work);
          work->scale(-0.5);
          num.add_block(dim, 0, dim, dim, work);
          num.add_block(0, dim, dim, dim, work);
          Matrix fss = *shalf_xh_ % num * *shalf_xh_;
          denom_xh_ = std::unique_ptr<double[]>(new double[2*dim]);
          fss.diagonalize(denom_xh_.get());
          *shalf_xh_ = fss % *shalf_xh_;
        }

        // xc/xx case
        {
          const size_t dim = nact*nact*nact;
          const size_t size = dim*dim;
          std::shared_ptr<RDM<3> > ovl(new RDM<3>(*rdm3source));
          ovl->scale(-1.0);
          for (int i5 = 0; i5 != nact; ++i5)
            for (int i4 = 0; i4 != nact; ++i4)
              for (int i2 = 0; i2 != nact; ++i2)
                for (int i3 = 0; i3 != nact; ++i3)
                  for (int i1 = 0; i1 != nact; ++i1)
                    for (int i0 = 0; i0 != nact; ++i0) {
                      if (i2 == i4)             ovl->element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2->element(i0, i1, i3, i5);
                      if (i2 == i3)             ovl->element(i0, i1, i3, i2, i4, i5) +=  2.0 * rdm2->element(i0, i1, i4, i5);
                      if (i1 == i4)             ovl->element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2->element(i3, i2, i0, i5);
                      if (i1 == i4 && i2 == i3) ovl->element(i0, i1, i3, i2, i4, i5) +=  2.0 * rdm1mat->element(i0, i5);
                      if (i1 == i3)             ovl->element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2->element(i0, i2, i4, i5);
                      if (i1 == i3 && i2 == i4) ovl->element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm1mat->element(i0, i5);
                    }
          shalf_xxh_ = std::shared_ptr<Matrix>(new Matrix(dim, dim));
          sort_indices<0,1,3,5,4,2,0,1,1,1>(ovl->data(), shalf_xxh_->data(), nact, nact, nact, nact, nact, nact);
          shalf_xxh_->inverse_half(1.0e-9);

          std::shared_ptr<RDM<4> > rdm4(new RDM<4>(*rdm4source));
          rdm4->scale(-1.0);
          for (int i4 = 0; i4 != nact; ++i4)
            for (int i3 = 0; i3 != nact; ++i3)
              for (int i7 = 0; i7 != nact; ++i7)
                for (int i6 = 0; i6 != nact; ++i6)
                  for (int i2 = 0; i2 != nact; ++i2)
                    for (int i5 = 0; i5 != nact; ++i5)
                      for (int i1 = 0; i1 != nact; ++i1)
                        for (int i0 = 0; i0 != nact; ++i0) {
                          if (i4 == i6)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i1, i5, i2, i3, i7);
                          if (i4 == i5)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i1, i3, i2, i6, i7);
                          if (i2 == i6)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i1, i3, i4, i5, i7);
                          if (i2 == i6 && i4 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i1, i3, i7);
                          if (i2 == i5)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm3source->element(i0, i1, i3, i4, i6, i7);
                          if (i2 == i5 && i4 == i6)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2->element(i0, i1, i3, i7);
                          if (i2 == i3)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i1, i5, i4, i6, i7);
                          if (i2 == i3 && i4 == i6)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i1, i5, i7);
                          if (i2 == i3 && i4 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2->element(i0, i1, i6, i7);
                          if (i1 == i6)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i3, i4, i5, i2, i0, i7);
                          if (i1 == i6 && i4 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i3, i2, i0, i7);
                          if (i1 == i6 && i2 == i3)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i5, i4, i0, i7);
                          if (i1 == i6 && i2 == i3 && i4 == i5) rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm1mat->element(i0, i7);
                          if (i1 == i5)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i2, i3, i4, i6, i7);
                          if (i1 == i5 && i4 == i6)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i2, i3, i7);
                          if (i2 == i3 && i1 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i4, i6, i7);
                          if (i2 == i3 && i1 == i5 && i4 == i6) rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm1mat->element(i0, i7);
                          if (i1 == i3)                         rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3source->element(i0, i4, i5, i2, i6, i7);
                          if (i1 == i3 && i4 == i6)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i5, i2, i0, i7);
                          if (i1 == i3 && i4 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i2, i6, i7);
                          if (i2 == i6 && i1 == i3)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i0, i4, i5, i7);
                          if (i2 == i6 && i1 == i3 && i4 == i5) rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm1mat->element(i0, i7);
                          if (i1 == i3 && i2 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2->element(i0, i4, i6, i7);
                          if (i4 == i6 && i2 == i5 && i1 == i3) rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm1mat->element(i0, i7);
                          if (i1 == i6 && i2 == i5)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2->element(i3, i4, i0, i7);
                          if (i1 == i5 && i2 == i6)             rdm4->element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2->element(i3, i4, i0, i7);
                        }
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, rdm4->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);
          Matrix fss(dim, dim);
          sort_indices<0,1,3,5,4,2,0,1,1,1>(work2->data(), fss.data(), nact, nact, nact, nact, nact, nact);
          fss = *shalf_xxh_ % fss * *shalf_xxh_;
          denom_xxh_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_xxh_.get());
          *shalf_xxh_ = fss % *shalf_xxh_; 
        }
        // xx/xa case
        {
          const size_t dim = nact*nact*nact;
          const size_t size = dim*dim;
          std::shared_ptr<RDM<3> > ovl(new RDM<3>(*rdm3source));
          for (int i5 = 0; i5 != nact; ++i5)
            for (int i0 = 0; i0 != nact; ++i0)
              for (int i4 = 0; i4 != nact; ++i4)
                for (int i3 = 0; i3 != nact; ++i3)
                  for (int i2 = 0; i2 != nact; ++i2)
                    for (int i1 = 0; i1 != nact; ++i1) {
                      if (i2 == i3) ovl->element(i1, i2, i3, i4, i0, i5) += 1.0 * rdm2->element(i1, i4, i0, i5);
                    }
          shalf_xhh_ = std::shared_ptr<Matrix>(new Matrix(dim, dim));
          sort_indices<4,0,1,5,3,2,0,1,1,1>(ovl->data(), shalf_xhh_->data(), nact, nact, nact, nact, nact, nact);
          shalf_xhh_->inverse_half(1.0e-9);

          std::shared_ptr<RDM<4> > rdm4(new RDM<4>(*rdm4source));
          for (int i4 = 0; i4 != nact; ++i4)
            for (int i3 = 0; i3 != nact; ++i3)
              for (int i7 = 0; i7 != nact; ++i7)
                for (int i0 = 0; i0 != nact; ++i0)
                  for (int i6 = 0; i6 != nact; ++i6)
                    for (int i5 = 0; i5 != nact; ++i5)
                      for (int i2 = 0; i2 != nact; ++i2)
                        for (int i1 = 0; i1 != nact; ++i1) {
                          if (i4 == i5)             rdm4->element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3source->element(i1, i2, i3, i6, i0, i7);
                          if (i2 == i3)             rdm4->element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3source->element(i1, i4, i5, i6, i0, i7);
                          if (i2 == i3 && i4 == i5) rdm4->element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm2->element(i1, i6, i0, i7);
                          if (i2 == i5)             rdm4->element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3source->element(i3, i4, i1, i6, i0, i7);
                        }
          std::shared_ptr<Matrix> work2(new Matrix(dim, dim));
          dgemv_("N", size, nact*nact, 1.0, rdm4->data(), size, fockact->data(), 1, 0.0, work2->data(), 1);
          Matrix fss(dim, dim);
          sort_indices<4,0,1,5,3,2,0,1,1,1>(work2->data(), fss.data(), nact, nact, nact, nact, nact, nact);
          fss = *shalf_xhh_ % fss * *shalf_xhh_;
          denom_xhh_ = std::unique_ptr<double[]>(new double[dim]);
          fss.diagonalize(denom_xhh_.get());
          *shalf_xhh_ = fss % *shalf_xhh_; 
        }
      }

      // set e0
      e0_ = compute_e0();
    }

    IndexRange& virt() { return virt_; }
    IndexRange& all() { return all_; }
    IndexRange& closed() { return closed_; }

    std::shared_ptr<const Reference>& ref() { return ref_; };

    double e0() const { return e0_; }

    virtual void solve() = 0;

};

}
}

#endif

