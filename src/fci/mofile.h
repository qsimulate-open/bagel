//
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


#ifndef __NEWINT_FCI_MOFILE_H
#define __NEWINT_FCI_MOFILE_H

#include <fstream>
#include <string>
#include <memory>
#include <cassert>
#include <src/wfn/reference.h>
#include <src/util/filename.h>
#include <src/scf/scf.h>
#include <src/scf/geometry.h>


class MOFile {

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

    std::unique_ptr<double[]> core_fock_;
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
    virtual double compute_mo1e(const int, const int) = 0;
    // this sets mo2e_1ext_ (half transformed DF integrals) and returns mo2e IN UNCOMPRESSED FORMAT
    virtual std::unique_ptr<double[]> compute_mo2e(const int, const int) = 0;

  public:
    MOFile(const std::shared_ptr<const Geometry>, const std::shared_ptr<const Reference>, const int, const int);
    ~MOFile();

    const std::shared_ptr<const Geometry> geom() const { return geom_; };

    const int sizeij() const { return sizeij_; };
    double mo1e(const size_t i) const { return mo1e_[i]; };
    double mo2e(const size_t i, const size_t j) const { return mo2e_[i+j*sizeij_]; };
    // strictly i <= j, k <= l
    double mo2e(const int i, const int j, const int k, const int l) const { return mo2e(address_(i,j), address_(k,l)); };
    double mo1e(const int i, const int j) const { return mo1e(address_(i,j)); };
    double* core_fock_ptr() { return core_fock_.get(); };
    double* mo1e_ptr() { return mo1e_.get(); };
    double* mo2e_ptr() { return mo2e_.get(); };
    const double* core_fock_ptr() const { return core_fock_.get(); };
    const double* mo1e_ptr() const { return mo1e_.get(); };
    const double* mo2e_ptr() const { return mo2e_.get(); };

    double core_energy() const { return core_energy_; };

    std::shared_ptr<DF_Half> mo2e_1ext() { return mo2e_1ext_; };
    std::shared_ptr<const DF_Half> mo2e_1ext() const { return mo2e_1ext_; };
    const double* const mo2e_1ext_ptr() const { return mo2e_1ext_->data(); };
    void update_1ext_ints(const std::vector<double>& coeff);

};

class Jop : public MOFile {
  protected:
    double compute_mo1e(const int, const int);
    std::unique_ptr<double[]> compute_mo2e(const int, const int);
  public:
    Jop(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const Reference> b, const int c, const int d) : MOFile(a,b,c,d) {
      core_energy_ = create_Jiiii(c, d);
    };
    ~Jop() {};
};

#endif
