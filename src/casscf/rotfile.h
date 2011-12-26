//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_CASSCF_ROTFILE_H
#define __NEWINT_CASSCF_ROTFILE_H

class RotFile {
  protected:
    double* data_;
    const int nclosed_;
    const int nact_;
    const int nvirt_;

  public:
    RotFile(const int iclos, const int iact, const int ivirt)
     : nclosed_(iclos), nact_(iact), nvirt_(ivirt) {
      data_ = new double[iclos*iact+iclos*ivirt+iact*ivirt];
    };

    ~RotFile() { delete[] data_; };

    // closed-active block. closed runs first
    double* ptr_ca() { return data_; };
    double& ptr_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; };
    // active-virtual block. virtual runs first
    double* ptr_va() { return data_ + nclosed_*nact_; };
    double& ptr_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; };
    // closed-virtual block. virtual runs first
    double* ptr_vc() { return data_ + (nclosed_+nvirt_)*nact_; };
    double& ptr_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; };

};

#endif
