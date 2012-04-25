//
// Author : Toru Shiozaki
// Date   : April 2012
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

  public:
    JKop(const std::shared_ptr<DensityFit> df, const std::shared_ptr<Coeff> c, const std::shared_ptr<Fock<1> > hcore,
         const std::shared_ptr<FCI> fci, const size_t nocc, const size_t nclosed, const size_t nact) : fci_(fci) {
      assert(nclosed+nact == nocc);
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
      std::fill(rdm2allk.get(), rdm2allk.get()+nocc*nocc*nocc*nocc, 0.0);
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
            for (int l = 0; l != nocc; ++l) { 
              rdm2allk[l+nocc*(k+nocc*(j+nocc*i))] = rdm2all[l+nocc*(j+nocc*(k+nocc*i))];
            }
          }
        }
      }
      // sort K 
      const size_t nbasis = df->nbasis();
      std::unique_ptr<double[]> buf(new double[nocc*nocc*nbasis*nbasis]);
      for (int i = 0; i != nocc; ++i) { 
        for (int j = 0; j != nocc; ++j) { 
          for (int k = 0; k != nbasis; ++k) {
            for (int l = 0; l != nbasis; ++l) {
              buf[l+nbasis*(k+nbasis*(j+nocc*i))] = data_[l+nbasis*(j+nocc*(k+nbasis*i))]; 
            }
          }
        }
      }
      // first K contribution to G operator
      dgemm_("N", "N", nbasis*nbasis, nocc*nocc, nocc*nocc, 2.0, buf, nbasis*nbasis, rdm2allk, nocc*nocc, 0.0, data_, nbasis*nbasis);
      // second J contribution to G operator
      dgemm_("N", "N", nbasis*nbasis, nocc*nocc, nocc*nocc, 1.0, jdata_, nbasis*nbasis, rdm2all, nocc*nocc, 0.0, data_, nbasis*nbasis);
      // last h contribution
      size_t icnt = 0lu;
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j, icnt += nbasis*nbasis) {
          daxpy_(nbasis*nbasis, rdm1_av->element(j,i), hcore->data(), 1, data_.get()+icnt, 1); 
        }
      }
    };
    ~JKop() {};

    std::shared_ptr<Matrix1e> contract(const std::shared_ptr<Matrix1e> in) {
#if 0
      std::shared_ptr<Matrix1e> out(fci_->geom());
      for (int i = 0; i != nocc; ++i) {
        for (int j = 0; j != nocc; ++j) {
          for (int k = 0; k != nbasis; ++k) {
            for (int l = 0; l != nbasis; ++l) {
              out->element(k, i) += in->element(l,j) * data_[l+nbasis*(k+nbasis*(j+nocc*i))]; 
            }
          }
        }
      }
      return out;
#endif
    };

};

#endif
