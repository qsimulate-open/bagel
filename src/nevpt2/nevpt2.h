//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_NEVPT2_NEVPT2_H
#define __SRC_NEVPT2_NEVPT2_H

#include <src/wfn/method.h>
#include <src/casscf/casscf.h>

namespace bagel {

class NEVPT2 : public Method {
  protected:
    std::shared_ptr<CASSCF> casscf_;
    int ncore_;
    int nclosed_;
    int nact_;
    int nvirt_;
    int istate_;
    double norm_thresh_;

    std::string abasis_;

    double energy_;

    // coefficient blocks
    std::shared_ptr<const Matrix> ccoeff_;
    std::shared_ptr<const Matrix> acoeff_;
    std::shared_ptr<const Matrix> vcoeff_;

    // density matrices to be used
    // particle RDMs
    std::shared_ptr<const Matrix> rdm1_;
    std::shared_ptr<const Matrix> rdm2_;
    std::shared_ptr<const Matrix> rdm3_;
    std::shared_ptr<const Matrix> rdm4_;
    // hole RDMs
    std::shared_ptr<const Matrix> hrdm1_;
    std::shared_ptr<const Matrix> hrdm2_;
    std::shared_ptr<const Matrix> hrdm3_;
    // <a+a b+b c+c..>
    std::shared_ptr<const Matrix> ardm2_;
    std::shared_ptr<const Matrix> ardm3_;
    std::shared_ptr<const Matrix> ardm4_;
    // <a+a bb+>
    std::shared_ptr<const Matrix> srdm2_;
    // <a+a bb+ c+c>
    std::shared_ptr<const Matrix> srdm3_;

    // integrals in physicists notation
    std::shared_ptr<const Matrix> ints2_;
    std::shared_ptr<const Matrix> fockact_;
    std::shared_ptr<const Matrix> fockact_c_;
    std::shared_ptr<const Matrix> fockact_h_;
    std::shared_ptr<const Matrix> fockact_p_;

    // K and K'mat
    std::shared_ptr<const Matrix> kmat_;
    std::shared_ptr<const Matrix> kmatp_;
    std::shared_ptr<const Matrix> kmat2_;
    std::shared_ptr<const Matrix> kmatp2_;

    // A
    std::shared_ptr<const Matrix> amat2_;
    std::shared_ptr<const Matrix> amat3_;
    std::shared_ptr<const Matrix> amat3t_;
    // B
    std::shared_ptr<const Matrix> bmat2_;
    std::shared_ptr<const Matrix> bmat2t_;
    // C
    std::shared_ptr<const Matrix> cmat2_;
    std::shared_ptr<const Matrix> cmat2t_;
    // D
    std::shared_ptr<const Matrix> dmat2_;
    std::shared_ptr<const Matrix> dmat1_;
    std::shared_ptr<const Matrix> dmat1t_;

    void compute_rdm();
    void compute_hrdm();
    void compute_asrdm();
    void compute_ints();
    void compute_kmat();
    void compute_abcd();

  public:
    NEVPT2(const std::shared_ptr<const PTree>, const std::shared_ptr<const Geometry>, const std::shared_ptr<const Reference> = nullptr);

    virtual void compute() override;
    virtual std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    double energy() const { return energy_; }
    int ncore() const { return ncore_; }
    std::string abasis() const { return abasis_; }
};

}

#endif
