//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mofile.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//



#ifndef __BAGEL_FCI_MOFILE_H
#define __BAGEL_FCI_MOFILE_H

#include <src/util/math/csymmatrix.h>
#include <src/scf/hf/rhf.h>

namespace bagel {

class MOFile {

  protected:
    int nocc_;

    bool hz_; // If true, do hz stuff. This may be revisited if more algorithms are implemented
    double core_energy_;

    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    size_t sizeij_;

    // core fock operator
    std::shared_ptr<const Matrix> core_fock_;
    // mo1e is a compressed symmetric matrix
    std::shared_ptr<CSymMatrix> mo1e_;
    // mo2e is a matrix of sizeij_*sizeij_
    std::shared_ptr<Matrix> mo2e_;

    std::shared_ptr<DFHalfDist> mo2e_1ext_;

    std::shared_ptr<const Matrix> coeff_;

    int address_(int i, int j) const { assert(i <= j); return i+((j*(j+1))>>1); }

    // creates integral files and returns the core energy.
    void init(const int nstart, const int nfence, const bool store);

    // this sets mo1e_
    virtual std::shared_ptr<const Matrix> compute_mo1e(const int, const int) = 0;
    // this sets mo2e_1ext_ (half transformed DF integrals) and returns mo2e IN UNCOMPRESSED FORMAT
    virtual std::shared_ptr<const Matrix> compute_mo2e(const int, const int) = 0;

    void compress_and_set(std::shared_ptr<const Matrix> buf1e, std::shared_ptr<const Matrix> buf2e);


  public:
    MOFile(const std::shared_ptr<const Reference>, const std::string method = std::string("KH"));
    MOFile(const std::shared_ptr<const Reference>, const std::shared_ptr<const Matrix>, const std::string method = std::string("KH"));
    // Shortcut used in MEH
    MOFile(const std::shared_ptr<CSymMatrix> mo1e, const std::shared_ptr<Matrix> mo2e) : nocc_(mo1e->nocc()), hz_(true), mo1e_(mo1e), mo2e_(mo2e) {}

    const std::shared_ptr<const Geometry> geom() const { return geom_; }

    int sizeij() const { return sizeij_; }
    std::shared_ptr<const CSymMatrix> mo1e() const { return mo1e_; }
    std::shared_ptr<const Matrix> mo2e() const { return mo2e_; }
    double& mo1e(const size_t i) { return mo1e_->data(i); }
    double& mo2e(const size_t i, const size_t j) { return mo2e_->element(i,j); }
    const double& mo1e(const size_t i) const { return mo1e_->data(i); }
    const double& mo2e(const size_t i, const size_t j) const { return mo2e_->element(i,j); }
    // input in (ij|kl), accesses the right way
    const double& mo2e(size_t i, size_t j, size_t k, size_t l) const { return (hz_ ? mo2e_hz(i, k, j, l) : mo2e_kh(i, j, k, l)); }
    // This is really ugly but will work until I can think of some elegant solution that keeps mo2e(i,j,k,l) inline but doesn't require more derived classes
    // strictly i <= j, k <= l
    double& mo2e_kh(const int i, const int j, const int k, const int l) { return mo2e(address_(i,j), address_(k,l)); }
    const double& mo2e_kh(const int i, const int j, const int k, const int l) const { return mo2e(address_(i,j), address_(k,l)); }
    // This is in <ij|kl> == (ik|jl) format
    double& mo2e_hz(const int i, const int j, const int k, const int l) { return mo2e_->element(i+nocc_*j, k+nocc_*l); }
    double& mo2e_hz(const int ij, const int kl) { return mo2e_->element(ij, kl); }
    double& mo1e(const int i, const int j) { return mo1e_->element(i,j); }
    const double& mo2e_hz(const int i, const int j, const int k, const int l) const { return mo2e_->element(i+nocc_*j, k+nocc_*l); }
    const double& mo2e_hz(const int ij, const int kl) const { return mo2e_->element(ij, kl); }
    const double& mo1e(const int i, const int j) const { return mo1e_->element(i,j); }
    double* mo1e_ptr() { return mo1e_->data(); }
    double* mo2e_ptr() { return mo2e_->data(); }
    const double* mo1e_ptr() const { return mo1e_->data(); }
    const double* mo2e_ptr() const { return mo2e_->data(); }

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }

    bool hz() const { return hz_; }
    int nocc() const { return nocc_; }

    std::shared_ptr<const Matrix> core_fock() const { return core_fock_; }
    double core_energy() const { return core_energy_; }

    std::shared_ptr<DFHalfDist> mo2e_1ext() { return mo2e_1ext_; }
    std::shared_ptr<const DFHalfDist> mo2e_1ext() const { return mo2e_1ext_; }
    void update_1ext_ints(const std::shared_ptr<const Matrix>& coeff);

};


class Jop : public MOFile {
  protected:
    std::shared_ptr<const Matrix> compute_mo1e(const int, const int) override;
    std::shared_ptr<const Matrix> compute_mo2e(const int, const int) override;
  public:
    Jop(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const Matrix> e, const bool store, const std::string f = "KH")
      : MOFile(b,e,f) { init(c, d, store); }
    Jop(const std::shared_ptr<CSymMatrix> mo1e, const std::shared_ptr<Matrix> mo2e) : MOFile(mo1e, mo2e) {}
};


class Htilde : public MOFile {
  protected:
    // temp storage
    std::shared_ptr<const Matrix> h1_tmp_;
    std::shared_ptr<const Matrix> h2_tmp_;

    std::shared_ptr<const Matrix> compute_mo1e(const int, const int) override { return h1_tmp_; }
    std::shared_ptr<const Matrix> compute_mo2e(const int, const int) override { return h2_tmp_; }

  public:
    Htilde(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const Matrix> h1, std::shared_ptr<const Matrix> h2)
      : MOFile(b), h1_tmp_(h1), h2_tmp_(h2) {
      assert(c == 0); // otherwise there will be an unnecessary computation
      init(c, d, false);
    }
};

}

#endif
