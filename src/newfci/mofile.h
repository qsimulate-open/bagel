
// Newint - Parallel electron correlation program.
// Filename: mofile.h
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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



#ifndef __NEWINT_NEWFCI_MOFILE_H
#define __NEWINT_NEWFCI_MOFILE_H

#define __NEWFCI_DEBUGGING

#include <fstream>
#include <string>
#include <memory>
#include <cassert>
#include <tuple>
#include <src/wfn/reference.h>
#include <src/util/filename.h>
#include <src/scf/scf.h>
#include <src/scf/geometry.h>


class NewMOFile {

  protected:
    int nocc_;
    int nbasis_;

    bool do_df_;
    double core_energy_;

    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const Reference> ref_;
    std::shared_ptr<std::fstream> file_;
    size_t sizeij_;
    long filesize_;
    std::string filename_;
    std::vector<std::shared_ptr<Shell> > basis_;
    std::vector<int> offset_;

    std::shared_ptr<Matrix1e> core_fock_;
    std::unique_ptr<double[]> mo1e_;
    std::unique_ptr<double[]> mo2e_;
    std::shared_ptr<DF_Half> mo2e_1ext_;

    size_t mo2e_1ext_size_;

    int address_(int i, int j) const {
      assert(i <= j);
      return i+((j*(j+1))>>1);
    };
    // creates integral files and returns the core energy.
    double create_Jiiii(const int, const int);

    // this sets mo1e_, core_fock_ and returns a core energy
    virtual std::tuple<std::unique_ptr<double[]>, double> compute_mo1e(const int, const int) {
      assert(false);
      std::unique_ptr<double[]> a; double b = 0.0; return std::make_tuple(std::move(a),b);
    };
    // this sets mo2e_1ext_ (half transformed DF integrals) and returns mo2e IN UNCOMPRESSED FORMAT
    virtual std::unique_ptr<double[]> compute_mo2e(const int, const int) { return std::unique_ptr<double[]>(); };

    void compress(std::unique_ptr<double[]>& buf1e, std::unique_ptr<double[]>& buf2e);

    std::shared_ptr<const Coeff> coeff_;

  public:
    NewMOFile(const std::shared_ptr<const Reference>, const int nstart, const int nfence);
    NewMOFile(const std::shared_ptr<const Reference>, const int nstart, const int nfence, const std::shared_ptr<const Coeff>);
    ~NewMOFile();

    const std::shared_ptr<const Geometry> geom() const { return geom_; };

    int sizeij() const { return sizeij_; };
    double mo1e(const size_t i) const { return mo1e_[i]; };
    double mo2e(const size_t i, const size_t j) const { return mo2e_[i+j*sizeij_]; };
    // strictly i <= j, k <= l
    #ifndef __NEWFCI_DEBUGGING
    double mo2e(const int i, const int j, const int k, const int l) const { return mo2e(address_(i,j), address_(k,l)); };
    #else
    double mo2e(const int i, const int j, const int k, const int l) const { return mo2e_[k + nocc_*j + nocc_*nocc_*l + nocc_*nocc_*nocc_*i]; };
    #endif
    double mo1e(const int i, const int j) const { return mo1e(address_(i,j)); };
    std::shared_ptr<const Matrix1e> core_fock() const { return core_fock_; };
    double* core_fock_ptr() { return core_fock_->data(); };
    double* mo1e_ptr() { return mo1e_.get(); };
    double* mo2e_ptr() { return mo2e_.get(); };
    const double* core_fock_ptr() const { return core_fock_->data(); };
    const double* mo1e_ptr() const { return mo1e_.get(); };
    const double* mo2e_ptr() const { return mo2e_.get(); };

    double core_energy() const { return core_energy_; };

    std::shared_ptr<DF_Half> mo2e_1ext() { return mo2e_1ext_; };
    std::shared_ptr<const DF_Half> mo2e_1ext() const { return mo2e_1ext_; };
    const double* mo2e_1ext_ptr() const { return mo2e_1ext_->data(); };
    void update_1ext_ints(const std::vector<double>& coeff);

};

class NewJop : public NewMOFile {
  protected:
    std::tuple<std::unique_ptr<double[]>, double> compute_mo1e(const int, const int);
    std::unique_ptr<double[]> compute_mo2e(const int, const int);
  public:
    NewJop(const std::shared_ptr<const Reference> b, const int c, const int d) : NewMOFile(b,c,d) { core_energy_ = create_Jiiii(c, d); };
    NewJop(const std::shared_ptr<const Reference> b, const int c, const int d, std::shared_ptr<const Coeff> e) : NewMOFile(b,c,d,e) { core_energy_ = create_Jiiii(c, d); };
    ~NewJop() {};
};

class NewHtilde : public NewMOFile {
  protected:
    // temp storage
    std::unique_ptr<double[]> h1_tmp_;
    std::unique_ptr<double[]> h2_tmp_;

    std::tuple<std::unique_ptr<double[]>, double> compute_mo1e(const int, const int) { return std::make_tuple(std::move(h1_tmp_), 0.0); }; 
    std::unique_ptr<double[]> compute_mo2e(const int, const int) { return std::move(h2_tmp_); };
  
  public:
    NewHtilde(const std::shared_ptr<const Reference> b, const int c, const int d, std::unique_ptr<double[]> h1, std::unique_ptr<double[]> h2)
      : NewMOFile(b,c,d), h1_tmp_(std::move(h1)), h2_tmp_(std::move(h2)) {
      core_energy_ = create_Jiiii(c, d);
    }
    ~NewHtilde() {};
};

#endif
