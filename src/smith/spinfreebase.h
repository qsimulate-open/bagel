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
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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
#include <src/smith/denom.h>
#include <src/math/matrix.h>
#include <src/fci/dvec.h>
#include <src/wfn/reference.h>
#include <src/prop/dipole.h>
#include <src/grad/cphf.h>
#include <chrono>

namespace bagel {
namespace SMITH {

template <typename T>
class SpinFreeMethod {
  protected:
    // deprecated
    IndexRange virt_;
    IndexRange active_;
    IndexRange closed_;
    IndexRange all_;
    IndexRange ci_;

    // new
    std::shared_ptr<const IndexRange> rvirt_;
    std::shared_ptr<const IndexRange> ractive_;
    std::shared_ptr<const IndexRange> rclosed_;
    std::shared_ptr<const IndexRange> rci_;

    std::shared_ptr<const Reference> ref_;

    std::shared_ptr<const Coeff> coeff_;
    double e0_;

    std::shared_ptr<Tensor<T>> v2_;
    std::shared_ptr<Tensor<T>> f1_;
    std::shared_ptr<Tensor<T>> h1_;
    std::shared_ptr<Tensor<T>> rdm1_;
    std::shared_ptr<Tensor<T>> rdm2_;
    std::shared_ptr<Tensor<T>> rdm3_;
    std::shared_ptr<Tensor<T>> rdm4_;

    // original determinants (for use in output)
    std::shared_ptr<const Determinants> det_;     

    // rdm ci derivatives
    std::shared_ptr<Tensor<T>> civec_;
    std::shared_ptr<Tensor<T>> rdm1deriv_;
    std::shared_ptr<Tensor<T>> rdm2deriv_;
    std::shared_ptr<Tensor<T>> rdm3deriv_;

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
      for (auto& i1 : active_) {
        for (auto& i0 : active_) {
          const size_t size = i0.size() * i1.size();
          std::unique_ptr<double[]> fdata = f1_->get_block(i0, i1);
          std::unique_ptr<double[]> rdata = rdm1_->get_block(i0, i1);
          sum += ddot_(size, fdata, 1, rdata, 1);
        }
      }
      std::cout << "    - Zeroth order energy: " << std::setw(20) << std::setprecision(10) << sum << std::endl;
      return sum;
    }

    std::shared_ptr<const Denom> denom_;
    std::shared_ptr<const Matrix> shalf_xhh() const { return denom_->shalf_xhh(); }
    std::shared_ptr<const Matrix> shalf_xxh() const { return denom_->shalf_xxh(); }
    std::shared_ptr<const Matrix> shalf_xh() const { return denom_->shalf_xh(); }
    std::shared_ptr<const Matrix> shalf_hh() const { return denom_->shalf_hh(); }
    std::shared_ptr<const Matrix> shalf_xx() const { return denom_->shalf_xx(); }
    std::shared_ptr<const Matrix> shalf_h() const { return denom_->shalf_h(); }
    std::shared_ptr<const Matrix> shalf_x() const { return denom_->shalf_x(); }
    const double& denom_xhh(const size_t i) const { return denom_->denom_xhh(i); }
    const double& denom_xxh(const size_t i) const { return denom_->denom_xxh(i); }
    const double& denom_xh(const size_t i) const { return denom_->denom_xh(i); }
    const double& denom_hh(const size_t i) const { return denom_->denom_hh(i); }
    const double& denom_xx(const size_t i) const { return denom_->denom_xx(i); }
    const double& denom_h(const size_t i) const { return denom_->denom_h(i); }
    const double& denom_x(const size_t i) const { return denom_->denom_x(i); }

    void update_amplitude(std::shared_ptr<Tensor<T>> t, const std::shared_ptr<Tensor<T>> r, const bool put = false) {

      // ranks of t and r are assumed to be the same

      // TODO should be parallelized
      for (auto& i3 : virt_) {
        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              if (!r->get_size(i0, i1, i2, i3)) continue;
              std::unique_ptr<double[]>       data0 = r->get_block(i0, i1, i2, i3);
              const std::unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);

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
                t->add_block(data0, i0, i1, i2, i3);
              } else {
                t->put_block(data0, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i2 : active_) {
        for (auto& i0 : active_) {
          // trans is the transformation matrix
          assert(shalf_xx());
          const int nact = ref_->nact();
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i0.size()*i2.size()*nact*nact]);
          for (int j2 = i2.offset(), k = 0; j2 != i2.offset()+i2.size(); ++j2)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
              std::copy_n(shalf_xx()->element_ptr(0,(j0-nclo)+(j2-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

          for (auto& i3 : virt_) {
            for (auto& i1 : virt_) {
              // if this block is not included in the current wave function, skip it
              if (!r->get_size(i0, i1, i2, i3)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(i0, i1, i2, i3);
              std::unique_ptr<double[]> data1(new double[r->get_size(i0, i1, i2, i3)]);
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
                    interm[iall] /= e0_ - (denom_xx(j02) + eig_[j3] + eig_[j1]);

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              dgemm_("T", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, transp, nact*nact, interm, nact*nact,
                                                                                    0.0, data0,  i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(data1, i0, i1, i2, i3);
              } else {
                t->put_block(data1, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i0 : active_) {
        // trans is the transformation matrix
        assert(shalf_x());
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        std::unique_ptr<double[]> transp(new double[i0.size()*nact]);
        for (int j0 = i0.offset(), k = 0; j0 != i0.offset()+i0.size(); ++j0, ++k)
          std::copy_n(shalf_x()->element_ptr(0,j0-nclo), nact, transp.get()+nact*k);

        for (auto& i3 : virt_) {
          for (auto& i2 : closed_) {
            for (auto& i1 : virt_) {
              if (!r->get_size(i2, i3, i0, i1)) continue;
              assert(r->get_size(i2, i1, i0, i3));
              std::unique_ptr<double[]>       data0 = r->get_block(i2, i3, i0, i1);
              const std::unique_ptr<double[]> data1 = r->get_block(i2, i1, i0, i3);
              std::unique_ptr<double[]> data2(new double[r->get_size(i2, i3, i0, i1)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
              sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());

              // move to orthogonal basis
              std::unique_ptr<double[]> interm(new double[i1.size()*i2.size()*i3.size()*nact]);
              dgemm_("N", "N", nact, i1.size()*i2.size()*i3.size(), i0.size(), 1.0, transp, nact, data2, i0.size(),
                                                                               0.0, interm, nact);

              size_t iall = 0;
              for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                  for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                    for (int j0 = 0; j0 != nact; ++j0, ++iall)
                      interm[iall] /= e0_ - (denom_x(j0) + eig_[j3] - eig_[j2] + eig_[j1]);

              // move back to non-orthogonal basis
              dgemm_("T", "N", i0.size(), i1.size()*i2.size()*i3.size(), nact, 1.0, transp, nact, interm, nact,
                                                                               0.0, data2,  i0.size());

              if (!put) {
                t->add_block(data2, i0, i1, i2, i3);
              } else {
                t->put_block(data2, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        // trans is the transformation matrix
        assert(shalf_h());
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        std::unique_ptr<double[]> transp(new double[i3.size()*nact]);
        for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3, ++k)
          std::copy_n(shalf_h()->element_ptr(0,j3-nclo), nact, transp.get()+nact*k);

        for (auto& i2 : closed_) {
          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              if (!r->get_size(i2, i3, i0, i1)) continue;
              assert(r->get_size(i0, i3, i2, i1));
              std::unique_ptr<double[]>       data0 = r->get_block(i2, i3, i0, i1);
              const std::unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);
              std::unique_ptr<double[]> data2(new double[r->get_size(i2, i3, i0, i1)]);
              sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
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
                      interm[iall] /= e0_ - (denom_h(j3) - eig_[j2] + eig_[j1] - eig_[j0]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i0.size()*i1.size()*i2.size(), i3.size(), nact, 1.0, interm, i0.size()*i1.size()*i2.size(), transp, nact,
                                                                               0.0, data2,  i0.size()*i1.size()*i2.size());

              if (!put) {
                t->add_block(data2, i0, i1, i2, i3);
              } else {
                t->put_block(data2, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          assert(shalf_hh());
          const int nact = ref_->nact();
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i1.size()*i3.size()*nact*nact]);
          for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++k)
              std::copy_n(shalf_hh()->element_ptr(0,(j1-nclo)+(j3-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

          for (auto& i2 : closed_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              if (!r->get_size(i0, i1, i2, i3)) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(i0, i1, i2, i3);
              std::unique_ptr<double[]> data1(new double[r->get_size(i0, i1, i2, i3)]);
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
                    interm[iall] /= e0_ - (denom_hh(j13) - eig_[j2] - eig_[j0]);

              // move back to non-orthogonal basis
              // factor of 0.5 due to the factor in the overlap
              dgemm_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, interm, i0.size()*i2.size(), transp, nact*nact,
                                                                                    0.0, data0,  i0.size()*i2.size());

              // sort back to the original order
              sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(data1, i0, i1, i2, i3);
              } else {
                t->put_block(data1, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          assert(shalf_xh());
          const int nact = ref_->nact();
          const int nclo = ref_->nclosed();
          std::unique_ptr<double[]> transp(new double[i2.size()*i3.size()*nact*nact*4]);
          for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++k) {
              std::copy_n(shalf_xh()->element_ptr(0,             (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*k);
              std::copy_n(shalf_xh()->element_ptr(0, nact*nact + (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*(k+i2.size()*i3.size()));
            }

          for (auto& i1 : virt_) {
            for (auto& i0 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->get_size(i2, i3, i0, i1);
              if (!blocksize) continue;
              assert(blocksize == r->get_size(i0, i3, i2, i1));
              std::unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
              std::unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);

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
                    interm[iall] /= e0_ - (denom_xh(j23) + eig_[j1] - eig_[j0]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i0.size()*i1.size(), i2.size()*i3.size()*2, nact*nact*2, 1.0, interm, i0.size()*i1.size(), transp, nact*nact*2,
                                                                                        0.0, data2,  i0.size()*i1.size());

              // sort back to the original order
              std::copy_n(data2.get(), blocksize, data0.get());
              sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksize, data1.get(), i0.size(), i1.size(), i2.size(), i3.size());
              if (!put) {
                t->add_block(data0, i0, i1, i2, i3);
                t->add_block(data1, i2, i1, i0, i3);
              } else {
                t->put_block(data0, i0, i1, i2, i3);
                t->put_block(data1, i2, i1, i0, i3);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i2 : active_) {
          for (auto& i0 : active_) {
            assert(shalf_xhh());
            const int nact = ref_->nact();
            const int nclo = ref_->nclosed();
            std::unique_ptr<double[]> transp(new double[i0.size()*i2.size()*i3.size()*nact*nact*nact]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                  std::copy_n(shalf_xhh()->element_ptr(0,j0-nclo+nact*(j2-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

            for (auto& i1 : virt_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->get_size(i2, i3, i0, i1);
              if (!blocksize) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
              std::unique_ptr<double[]> data1(new double[blocksize]);
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
                  interm[iall] /= e0_ - (denom_xhh(j123) + eig_[j1]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i1.size(), i0.size()*i2.size()*i3.size(), nact*nact*nact, 1.0, interm, i1.size(), transp, nact*nact*nact,
                                                                                         0.0, data0,  i1.size());

              // sort back to the original order
              sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i1.size(), i0.size(), i2.size(), i3.size());
              if (!put) {
                t->add_block(data1, i0, i1, i2, i3);
              } else {
                t->put_block(data1, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      for (auto& i3 : active_) {
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            assert(shalf_xxh());
            const int nact = ref_->nact();
            const int nclo = ref_->nclosed();
            std::unique_ptr<double[]> transp(new double[i0.size()*i1.size()*i3.size()*nact*nact*nact]);
            for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
                  std::copy_n(shalf_xxh()->element_ptr(0,j0-nclo+nact*(j1-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

            for (auto& i2 : closed_) {
              // if this block is not included in the current wave function, skip it
              const size_t blocksize = r->get_size(i2, i3, i0, i1);
              if (!blocksize) continue;
              // data0 is the source area
              std::unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
              std::unique_ptr<double[]> data1(new double[blocksize]);
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
                  interm[iall] /= e0_ - (denom_xxh(j013) - eig_[j2]);

              // move back to non-orthogonal basis
              dgemm_("N", "N", i2.size(), i0.size()*i1.size()*i3.size(), nact*nact*nact, 1.0, interm, i2.size(), transp, nact*nact*nact,
                                                                                         0.0, data0,  i2.size());

              // sort back to the original order
              sort_indices<1,2,0,3,0,1,1,1>(data0, data1, i2.size(), i0.size(), i1.size(), i3.size());
              if (!put) {
                t->add_block(data1, i0, i1, i2, i3);
              } else {
                t->put_block(data1, i0, i1, i2, i3);
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

      std::shared_ptr<const Dvec> dci0 = r->civectors();
      // TODO this should be updated when reference has civec
      std::shared_ptr<const Civec> bagel_civec = dci0->data(0);
      det_ = bagel_civec->det();

      // length of the ci expansion
      const size_t ci_size = r->civectors()->data(0)->size();
      ci_ = IndexRange(ci_size, max);

      rclosed_ = std::make_shared<const IndexRange>(c);
      ractive_ = std::make_shared<const IndexRange>(act);
      rvirt_   = std::make_shared<const IndexRange>(v);
      rci_     = std::make_shared<const IndexRange>(ci_);

      // f1 tensor.
      {
        std::vector<IndexRange> o = {all_, all_};
        MOFock<T> fock(ref_, o);
        f1_ = fock.tensor();
        h1_ = fock.hcore();
        // canonical orbitals within closed and virtual subspaces
        coeff_ = fock.coeff();
      }

      // v2 tensor.
      {
        IndexRange occ(closed_);
        occ.merge(active_);
        IndexRange virt(active_);
        virt.merge(virt_);

        std::vector<IndexRange> o = {occ, virt, occ, virt};
        K2ext<T> v2k(ref_, coeff_, o);
        v2_ = v2k.tensor();
      }

      // make a ci tensor.
      {
        std::vector<IndexRange> o = {ci_}; 
        // TODO fix later when referece has civec
        Ci<T> dci(ref_, o, bagel_civec);
        civec_ = dci.tensor();
      }

      // rdm ci derivatives. 
      {
        std::shared_ptr<const Dvec> rdm1d = r->rdm1deriv();

        std::vector<IndexRange> o = {ci_, active_, active_};
        rdm1deriv_ = std::make_shared<Tensor<T>>(o, false);
        const int nclo = ref_->nclosed();
        for (auto& i0 : active_) {
          for (auto& i1 : active_) {
            for (auto& ci0 : ci_) {
              const size_t size = i0.size() * i1.size() * ci0.size();
              std::unique_ptr<double[]> data(new double[size]);
              int iall = 0;
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is creation
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
                  for (int j2 = ci0.offset(); j2 != ci0.offset()+ci0.size(); ++j2, ++iall)
                    // Dvec - first index is annihilation, second is creation (see const_phis_ in fci/determinants.h and knowles_compute.cc)
                    data[iall] = rdm1d->data((j1-nclo)+r->nact()*(j0-nclo))->data(j2); 
              rdm1deriv_->put_block(data, ci0, i1, i0);
            }
          }
        }
      }

      {
        std::shared_ptr<const Dvec> rdm2d = r->rdm2deriv();

        std::vector<IndexRange> o = {ci_, active_, active_, active_, active_};
        rdm2deriv_ = std::make_shared<Tensor<T>>(o, false);
        const int nclo = ref_->nclosed();
        for (auto& i0 : active_) {
          for (auto& i1 : active_) {
            for (auto& i2 : active_) {
              for (auto& i3 : active_) {
                for (auto& ci0 : ci_) {
                  const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * ci0.size();
                  std::unique_ptr<double[]> data(new double[size]);
                  int iall = 0;
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is creation
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
                      for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) // this is creation
                        for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3) // this is annihilation
                          for (int j4 = ci0.offset(); j4 != ci0.offset()+ci0.size(); ++j4, ++iall)
                            data[iall] = rdm2d->data((j3-nclo)+r->nact()*((j2-nclo)+r->nact()*((j1-nclo)+r->nact()*(j0-nclo))))->data(j4); 
                  rdm2deriv_->put_block(data, ci0, i3, i2, i1, i0);
                }
              }
            }
          }
        }
      }

      {
        std::shared_ptr<const Dvec> rdm3d = r->rdm3deriv();

        std::vector<IndexRange> o = {ci_, active_, active_, active_, active_, active_, active_};
        rdm3deriv_ = std::make_shared<Tensor<T>>(o, false);
        const int nclo = ref_->nclosed();
        for (auto& i0 : active_) {
          for (auto& i1 : active_) {
            for (auto& i2 : active_) {
              for (auto& i3 : active_) {
                for (auto& i4 : active_) {
                  for (auto& i5 : active_) {
                    for (auto& ci0 : ci_) {
                      const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size() * ci0.size();
                      std::unique_ptr<double[]> data(new double[size]);
                      int iall = 0;
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is  creation
                        for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
                          for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) // this is creation
                            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3) // this is annihilation
                              for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4) // this is creation
                                for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5) // this is annhilation
                                  for (int j6 = ci0.offset(); j6 != ci0.offset()+ci0.size(); ++j6, ++iall)
                                    data[iall] = rdm3d->data((j5-nclo)+r->nact()*((j4-nclo)+r->nact()*((j3-nclo)+r->nact()*((j2-nclo)+r->nact()*((j1-nclo)+r->nact()*((j0-nclo)))))))->data(j6); 
                      rdm3deriv_->put_block(data, ci0, i5, i4, i3, i2, i1, i0);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // rdms.
      if (!ref_->rdm1().empty()) {
        std::vector<IndexRange> o = {active_, active_};
        rdm1_ = std::make_shared<Tensor<T>>(o, false);
        const int nclo = ref_->nclosed();
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            const size_t size = i0.size() * i1.size();
            std::unique_ptr<double[]> data(new double[size]);
            int iall = 0;
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                // TODO for the time being we hardwire "0" here (but this should be fixed)
                data[iall] = ref_->rdm1(0)->element(j0-nclo, j1-nclo);
            rdm1_->put_block(data, i0, i1);
          }
        }
      }
      if (!ref_->rdm2().empty()) {
        std::vector<IndexRange> o = {active_, active_, active_, active_};
        rdm2_ = std::make_shared<Tensor<T>>(o, false);
        const int nclo = ref_->nclosed();
        for (auto& i3 : active_) {
          for (auto& i2 : active_) {
            for (auto& i1 : active_) {
              for (auto& i0 : active_) {
                const size_t size = i0.size() * i1.size() * i2.size() * i3.size();
                std::unique_ptr<double[]> data(new double[size]);
                int iall = 0;
                for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                      for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                        // TODO for the time being we hardwire "0" here (but this should be fixed)
                        data[iall] = ref_->rdm2(0)->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo);
                 rdm2_->put_block(data, i0, i1, i2, i3);
              }
            }
          }
        }
      }
      // TODO generic function??
      if (!ref_->rdm1().empty() && !ref_->rdm2().empty()) {
        {
          std::vector<IndexRange> o = {active_, active_, active_, active_, active_, active_};
          rdm3_ = std::make_shared<Tensor<T>>(o, false);
          std::vector<IndexRange> p = {active_, active_, active_, active_, active_, active_, active_, active_};
          rdm4_ = std::make_shared<Tensor<T>>(p, false);
        }

        // TODO for the time being we hardwire "0" here (but this should be fixed)
        std::shared_ptr<RDM<3>> rdm3;
        std::shared_ptr<RDM<4>> rdm4;
        std::tie(rdm3, rdm4) = ref_->compute_rdm34(0);

        const int nclo = ref_->nclosed();
        for (auto& i5 : active_)
          for (auto& i4 : active_)
            for (auto& i3 : active_)
              for (auto& i2 : active_)
                for (auto& i1 : active_)
                  for (auto& i0 : active_) {
                    const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size();
                    std::unique_ptr<double[]> data(new double[size]);
                    int iall = 0;
                    for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                      for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                        for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                          for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                                data[iall] = rdm3->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo);
                    rdm3_->put_block(data, i0, i1, i2, i3, i4, i5);
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
                                        data[iall] = rdm4->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo, j6-nclo, j7-nclo);
                    rdm4_->put_block(data, i0, i1, i2, i3, i4, i5, i6, i7);
                  }


        const int nact = ref_->nact();
        auto fockact = std::make_shared<Matrix>(nact, nact);
        for (auto& i1 : active_)
          for (auto& i0 : active_)
            fockact->copy_block(i0.offset()-nclo, i1.offset()-nclo, i0.size(), i1.size(), this->f1_->get_block(i0, i1));

        // TODO hardwired 0
        auto rdm1 = std::make_shared<RDM<1>>(*ref_->rdm1(0));
        auto rdm2 = std::make_shared<RDM<2>>(*ref_->rdm2(0));

        // construct denominator
        denom_ = std::make_shared<const Denom>(*rdm1, *rdm2, *rdm3, *rdm4, *fockact);

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

    std::shared_ptr<const Civec> civec() const {
      return civec_->civec(det_);
    }

    Dipole dipole(std::shared_ptr<const Matrix> dm1, double correction) const {
      const size_t nclo = ref_->nclosed();
      const size_t nact = ref_->nact();

      // compute unrelaxed dipole moment
      // total density matrix
      auto dtot = std::make_shared<Matrix>(*dm1);

      // add correction to active space
      for (int i = nclo; i != nclo+nact; ++i) dtot->element(i,i) -=  correction*2.0;
      dtot->print("dm1 post correction", 20);

      for (int i = 0; i != nclo; ++i) dtot->element(i,i) += 2.0;
      // add to active space
      dtot->add_block(1.0, nclo, nclo, nact, nact, ref_->rdm1(0)->data());
      // convert to ao basis
      auto dtotao = std::make_shared<Matrix>(*coeff_ * *dtot ^ *coeff_);
      Dipole dipole(ref_->geom(), dtotao);
      return dipole;

    }


};

}
}

#endif

