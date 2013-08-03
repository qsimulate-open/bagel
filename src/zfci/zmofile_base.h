
// BAGEL - Parallel electron correlation program.
// Filename: zmofile_base.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell  <caldwell@u.northwestern.edu>
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



#ifndef __BAGEL_ZFCI_ZMOFILE_BASE_H
#define __BAGEL_ZFCI_ZMOFILE_BASE_H

#include <src/wfn/reference.h>
#include <src/math/zmatrix.h>

namespace bagel {

class ZMOFile_Base {

  protected:
    int nocc_;
    int nbasis_;

    bool do_df_;
    bool hz_; // If true, do hz stuff. This may be revisited if more algorithms are implemented
    double core_energy_;

    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const Reference> ref_;
    size_t sizeij_;
#if 0
    long filesize_;
    std::string filename_;
    std::vector<std::shared_ptr<Shell>> basis_;
    std::vector<int> offset_;
#endif
    std::unique_ptr<std::complex<double>[]> mo1e_;
    std::unique_ptr<std::complex<double>[]> mo2e_;
    std::shared_ptr<DFHalfDist> mo2e_1ext_;

    size_t mo2e_1ext_size_;

    // creates integral files and returns the core energy.
    virtual double create_Jiiii(const int, const int) = 0;
    // this sets mo1e_, core_fock_ and returns a core energy
    virtual std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) = 0;
    // this sets mo2e_1ext_ (half transformed DF integrals) and returns mo2e IN UNCOMPRESSED FORMAT
    virtual std::unique_ptr<std::complex<double>[]> compute_mo2e(const int, const int) = 0;
    virtual void compress(std::shared_ptr<const ZMatrix> buf1e, std::unique_ptr<std::complex<double>[]>& buf2e) = 0;
    std::shared_ptr<const Coeff> coeff_;
  public:
    ZMOFile_Base(const std::shared_ptr<const Reference> ref, const std::string method = std::string("KH")) :
      geom_(ref->geom()), ref_(ref), coeff_(ref_->coeff()) { }
    ZMOFile_Base(const std::shared_ptr<const Reference> ref, const std::shared_ptr<const Coeff> c, const std::string method = std::string("KH")) :
      hz_(false), geom_(ref->geom()), ref_(ref), coeff_(c) { }

    const std::shared_ptr<const Geometry> geom() const { return geom_; };

    int sizeij() const { return sizeij_; };

    void set_moints(std::shared_ptr<ZMatrix> mo1e, std::unique_ptr<std::complex<double>[]>& mo2e) {
      std::copy_n(mo1e->data(), nocc_*nocc_, mo1e_.get());
      mo2e_ = std::move(mo2e);
    }

    std::complex<double> mo1e(const size_t i) const { return mo1e_[i]; };
    std::complex<double> mo2e(const size_t i, const size_t j) const { return mo2e_[i+j*sizeij_]; };
    std::complex<double> mo1e(const size_t i, const size_t j) const { return mo1e_[i+j*nocc_]; };
    std::complex<double> mo2e_kh(const int i, const int j, const int k, const int l) const { return mo2e_[i+nocc_*(j+nocc_*(k+nocc_*l))]; };
    // This is in <ij|kl> == (ik|jl) format
#if 0
    std::complex<double> mo2e_hz(const int i, const int j, const int k, const int l) const { return mo2e_[i+nocc_*(j+nocc_*(k+nocc_*l))]; };
#endif
    std::complex<double>* mo1e_ptr() { return mo1e_.get(); };
    std::complex<double>* mo2e_ptr() { return mo2e_.get(); };
    const std::complex<double>* mo1e_ptr() const { return mo1e_.get(); };
    const std::complex<double>* mo2e_ptr() const { return mo2e_.get(); };

    double core_energy() const { return core_energy_; };
#if 0
    std::shared_ptr<DFHalfDist> mo2e_1ext() { return mo2e_1ext_; };
    std::shared_ptr<const DFHalfDist> mo2e_1ext() const { return mo2e_1ext_; };
    void update_1ext_ints(const std::shared_ptr<const Matrix>& coeff);
#endif

    void print(const int n = 10) const {
      for (int i = 0; i != std::min(nocc_,n); ++i)
        for (int j = 0; j != std::min(nocc_,n); ++j)
          for (int k = 0; k != std::min(nocc_,n); ++k)
            for (int l = 0; l != std::min(nocc_,n); ++l)
              std::cout << i << " " << j << " " << k << " " << l << " " << std::setprecision(8) << mo2e_kh(i,j,k,l) << std::endl;
    }
};
# if 0
class ZJop : public ZMOFile {
  protected:
    std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) override;
    std::unique_ptr<std::complex<double>[]> compute_mo2e(const int, const int) override;
  public:
    ZJop(const std::shared_ptr<const Reference> b, const int c, const int d, const std::string e = std::string("KH"))
      : ZMOFile(b,e) { core_energy_ = create_Jiiii(c, d); assert(false); }
    ZJop(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const Coeff> e, const std::string f = std::string("KH"))
      : ZMOFile(b,e,f) { core_energy_ = create_Jiiii(c, d); }
};
class ZHtilde : public ZMOFile {
  protected:
    // temp storage
    std::shared_ptr<const ZMatrix> h1_tmp_;
    std::unique_ptr<std::complex<double>[]> h2_tmp_;

    std::tuple<std::shared_ptr<const ZMatrix>, double> compute_mo1e(const int, const int) override { return std::make_tuple(h1_tmp_, 0.0); };
    std::unique_ptr<std::complex<double>[]> compute_mo2e(const int, const int) override { return std::move(h2_tmp_); };

  public:
    ZHtilde(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const ZMatrix> h1, std::unique_ptr<std::complex<double>[]> h2)
      : ZMOFile(b), h1_tmp_(h1), h2_tmp_(std::move(h2)) {
      core_energy_ = create_Jiiii(c, d);
    }
};
#endif
}

#endif
