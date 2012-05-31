//
// Newint - Parallel electron correlation program.
// Filename: jkop.h
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


#ifndef __SRC_CASSCF_JKOP_H
#define __SRC_CASSCF_JKOP_H

// implements external operators.
#include <src/fci/fci.h>
#include <algorithm>
#include <src/util/f77.h>

class JKop {
  protected:
    // CAUTION packing is different between J and K ops
    std::unique_ptr<double[]> jdata_;
    std::unique_ptr<double[]> data_;
    const std::shared_ptr<FCI> fci_;
    const size_t nocc_;
    const size_t nclosed_;
    const size_t nbasis_;

  public:
    JKop(const std::shared_ptr<const DensityFit> df, const std::shared_ptr<const Coeff> c, const std::shared_ptr<Fock<1> > hcore,
         const std::shared_ptr<FCI> fci, const size_t nocc, const size_t nclosed, const size_t nact)
    : fci_(fci), nocc_(nocc), nclosed_(nclosed), nbasis_(df->nbasis0()) {

      assert(nclosed+nact == nocc);
      assert(df->nbasis0() == df->nbasis1());
      // K operator
      std::shared_ptr<DF_Half> half = df->compute_half_transform(c->data(), nocc)->apply_J();
      data_ = std::move(half->form_4index());

      // J operator
      std::shared_ptr<DF_Full> full = half->compute_second_transform(c->data(), nocc);
      jdata_ = std::move(full->form_4index(df));

      // contruct 2RDM
      std::shared_ptr<RDM<1> > rdm1_av = fci->rdm1_av();
      std::shared_ptr<RDM<2> > rdm = fci->rdm2_av();
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
            dcopy_(nact, rdm->data()+nact*(k+nact*(j+nact*i)),1, rdm2all.get()+nclosed+nocc*(k+nclosed+nocc*(j+nclosed+nocc*(i+nclosed))),1);
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
            dcopy_(nocc, rdm2all.get()+nocc*(j+nocc*(k+nocc*i)), 1, rdm2allk.get()+nocc*(k+nocc*(j+nocc*i)), 1);
          }
        }
      }
      // sort K 
      std::unique_ptr<double[]> buf(new double[nocc*nocc*nbasis_*nbasis_]);
      for (int i = 0; i != nocc; ++i) { 
        for (int j = 0; j != nocc; ++j) { 
          for (int k = 0; k != nbasis_; ++k) {
            for (int l = 0; l != nbasis_; ++l) {
              buf[l+nbasis_*(k+nbasis_*(j+nocc*i))] = data_[j+nocc*(l+nbasis_*(i+nocc*k))]; 
            }
          }
        }
      }
      // first K contribution to G operator
      dgemm_("N", "N", nbasis_*nbasis_, nocc*nocc, nocc*nocc, 2.0, buf, nbasis_*nbasis_, rdm2allk, nocc*nocc, 0.0, data_, nbasis_*nbasis_);
      // second J contribution to G operator
      dgemm_("N", "N", nbasis_*nbasis_, nocc*nocc, nocc*nocc, 1.0, jdata_, nbasis_*nbasis_, rdm2all, nocc*nocc, 1.0, data_, nbasis_*nbasis_);
      // last h contribution
      size_t icnt = 0lu;
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j, icnt += nbasis_*nbasis_) {
          if (i >= nclosed && j >= nclosed) { 
            daxpy_(nbasis_*nbasis_, rdm1_av->element(j-nclosed,i-nclosed), hcore->data(), 1, data_.get()+icnt, 1); 
          } else if (i == j) {
            daxpy_(nbasis_*nbasis_, 2.0, hcore->data(), 1, data_.get()+icnt, 1); 
          } 
        }
      }
    };
    ~JKop() {};

    std::shared_ptr<Matrix1e> contract(const std::shared_ptr<Matrix1e> in) {
      std::shared_ptr<Matrix1e> out(new Matrix1e(fci_->geom()));
      const int nocc = nocc_;
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j) {
          for (int k = 0; k != nbasis_; ++k) {
            out->element(k, i) += 2.0*ddot_(nbasis_, in->element_ptr(0,j), 1, data_.get()+nbasis_*(k+nbasis_*(j+nocc*i)), 1);
          }
        }
      }
      return out;
    };

    std::shared_ptr<Matrix1e> denom() const {
      std::shared_ptr<Matrix1e> out(new Matrix1e(fci_->geom()));
      const size_t nbasis = nbasis_; 
      const int nocc = nocc_;
      const int nclosed = nclosed_;
      // virtual-occ part
      out->fill(1.0e1);
      for (int i = 0; i != nocc; ++i) {
        for (int j = nocc; j != nbasis; ++j) {
          out->element(j, i) = out->element(i, j) = 2.0*data_[j+nbasis*(j+nbasis*(i+nocc*i))];
        }
      }
      // occ-occ part
      for (int i = 0; i != nclosed; ++i) {
        for (int j = nclosed; j != nocc; ++j) {
          out->element(j, i) = out->element(i, j) = 2.0*data_[j+nbasis*(j+nbasis*(i+nocc*i))] + 2.0*data_[i+nbasis*(i+nbasis*(j+nocc*j))]
                                                    - 4.0*data_[j+nbasis*(i+nbasis*(j+nocc*i))];
        }
      }
      return out;
    };

};

#endif
