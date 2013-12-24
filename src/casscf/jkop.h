//
// BAGEL - Parallel electron correlation program.
// Filename: jkop.h
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


#ifndef __SRC_CASSCF_JKOP_H
#define __SRC_CASSCF_JKOP_H

// implements external operators.
#include <src/fci/fci.h>

namespace bagel {

class JKop {
  protected:
    // CAUTION packing is different between J and K ops
    std::shared_ptr<Matrix> jdata_;
    std::shared_ptr<Matrix> data_;
    const std::shared_ptr<FCI> fci_;
    const std::shared_ptr<const Coeff> coeff_;
    const size_t nocc_;
    const size_t nclosed_;
    const size_t nbasis_;

  public:
    JKop(const std::shared_ptr<const DFDist> df, const std::shared_ptr<const Coeff> c, const std::shared_ptr<const Matrix> hcore,
         const std::shared_ptr<FCI> fci, const size_t nocc, const size_t nclosed, const size_t nact)
    : fci_(fci), coeff_(c), nocc_(nocc), nclosed_(nclosed), nbasis_(df->nbasis0()) {

      std::shared_ptr<const Matrix> ocoeff = coeff_->slice(0, nocc);
      assert(nclosed+nact == nocc);
      assert(df->nbasis0() == df->nbasis1());
      // K operator // (ai|ai)
      std::shared_ptr<DFHalfDist> half = df->compute_half_transform(ocoeff)->apply_J();
      data_ = half->form_4index(half, 1.0);

      // J operator // (aa|ii)
      std::shared_ptr<DFFullDist> full = half->compute_second_transform(ocoeff);
      jdata_ = full->form_4index(df, 1.0);

      // contruct 2RDM
      std::shared_ptr<RDM<1>> rdm1_av = fci->rdm1_av();
      std::shared_ptr<RDM<2>> rdm = fci->rdm2_av();
      std::unique_ptr<double[]> rdm2all(new double[nocc*nocc*nocc*nocc]);
      std::unique_ptr<double[]> rdm2allk(new double[nocc*nocc*nocc*nocc]);
      std::fill(rdm2all.get(), rdm2all.get()+nocc*nocc*nocc*nocc, 0.0);
      // closed-closed
      for (int i = 0; i != nclosed; ++i) {
        for (int j = 0; j != nclosed; ++j) {
          rdm2all[j+nocc*(j+nocc*(i+nocc*i))] += 4.0;
          rdm2all[j+nocc*(i+nocc*(i+nocc*j))] -= 2.0;
        }
      }
      // active-active
      for (int i = 0; i != nact; ++i) {
        for (int j = 0; j != nact; ++j) {
          for (int k = 0; k != nact; ++k) {
            std::copy_n(rdm->data()+nact*(k+nact*(j+nact*i)), nact, rdm2all.get()+nclosed+nocc*(k+nclosed+nocc*(j+nclosed+nocc*(i+nclosed))));
          }
        }
      }
      // active-closed
      for (int i = 0; i != nclosed; ++i) {
        for (int j = 0; j != nact; ++j) {
          for (int k = 0; k != nact; ++k) {
            rdm2all[i+nocc*(i+nocc*(k+nclosed+nocc*(j+nclosed)))] += rdm1_av->element(k,j) * 2.0;
            rdm2all[k+nclosed+nocc*(j+nclosed+nocc*(i+nocc*(i)))] += rdm1_av->element(k,j) * 2.0;
            rdm2all[i+nocc*(k+nclosed+nocc*(j+nclosed+nocc*(i)))] -= rdm1_av->element(k,j);
            rdm2all[k+nclosed+nocc*(i+nocc*(i+nocc*(j+nclosed)))] -= rdm1_av->element(k,j);
          }
        }
      }
      // sort RDM for K
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j) {
          for (int k = 0; k != nocc; ++k) {
            std::copy_n(rdm2all.get()+nocc*(j+nocc*(k+nocc*i)), nocc, rdm2allk.get()+nocc*(k+nocc*(j+nocc*i)));
          }
        }
      }
      // sort K
      // after here buf contains K as (aa|ii)
      std::unique_ptr<double[]> buf(new double[nocc*nocc*nbasis_*nbasis_]);
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j) {
          for (int k = 0; k != nbasis_; ++k) {
            for (int l = 0; l != nbasis_; ++l) {
              buf[l+nbasis_*(k+nbasis_*(j+nocc*i))] = data_->element(j+nocc*l,i+nocc*k);
            }
          }
        }
      }
      // first K contribution to G operator
      dgemm_("N", "N", nbasis_*nbasis_, nocc*nocc, nocc*nocc, 2.0, buf.get(), nbasis_*nbasis_, rdm2allk.get(), nocc*nocc, 0.0, data_->data(), nbasis_*nbasis_);
      // second J contribution to G operator
      dgemm_("N", "N", nbasis_*nbasis_, nocc*nocc, nocc*nocc, 1.0, jdata_->data(), nbasis_*nbasis_, rdm2all.get(), nocc*nocc, 1.0, data_->data(), nbasis_*nbasis_);
      // last h contribution
      size_t icnt = 0lu;
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j, icnt += nbasis_*nbasis_) {
          if (i >= nclosed && j >= nclosed) {
            blas::ax_plus_y_n(rdm1_av->element(j-nclosed,i-nclosed), hcore->data(), nbasis_*nbasis_, data_->data()+icnt);
          } else if (i == j) {
            blas::ax_plus_y_n(2.0, hcore->data(), nbasis_*nbasis_, data_->data()+icnt);
          }
        }
      }
    }
    ~JKop() {}

    std::shared_ptr<Matrix> contract(const std::shared_ptr<Matrix> in) {
      auto out = std::make_shared<Matrix>(nbasis_, nbasis_);
      const int nocc = nocc_;
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j) {
          for (int k = 0; k != nbasis_; ++k) {
            out->element(k, i) += 2.0*ddot_(nbasis_, in->element_ptr(0,j), 1, data_->data()+nbasis_*(k+nbasis_*(j+nocc*i)), 1);
          }
        }
      }
      return out;
    }

    std::shared_ptr<Matrix> denom() const {
      const size_t nbasis = nbasis_;
      auto out = std::make_shared<Matrix>(nbasis, nbasis);
      const int nocc = nocc_;
      const int nclosed = nclosed_;
      auto tmp = std::make_shared<Matrix>(nbasis, nocc);
      // TODO this is an awful code.
      // first transform to MO
      {
        std::unique_ptr<double[]> buf(new double[nocc*nocc*nbasis*nbasis]);
        dgemm_("T", "N", nbasis, nocc*nocc*nbasis, nbasis, 1.0, coeff_->data(), nbasis, data_->data(), nbasis, 0.0, buf.get(), nbasis);
        for (int i = 0; i != nocc*nocc; ++i) {
          dgemm_("N", "N", nbasis, nbasis, nbasis, 1.0, buf.get()+i*nbasis*nbasis, nbasis, coeff_->data(), nbasis,
                                                   0.0, data_->data()+i*nbasis*nbasis, nbasis);
        }
      }

      // virtual-occ part
      out->fill(1.0e100);
      for (int i = 0; i != nocc; ++i) {
        for (int j = nocc; j != nbasis; ++j) {
          out->element(j, i) = out->element(i, j) = 2.0*jdata_->element(j+nbasis*j, i+nocc*i); // FIXME -- is this right?? looks as if it accesses to Jop
        }
      }
      // occ-occ part
      for (int i = 0; i != nclosed; ++i) {
        for (int j = nclosed; j != nocc; ++j) {
          out->element(j, i) = out->element(i, j) = (2.0*data_->element(j+nbasis*j, i+nocc*i) - 2.0*data_->element(i+nbasis*i, j+nocc*j));
        }
      }
      return out;
    }

};

}

#endif
