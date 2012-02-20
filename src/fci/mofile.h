//
// Newint - Parallel electron correlation program.
// Filename: mofile.h
// Copyright (C) 2011 Toru Shiozaki
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

    const std::shared_ptr<Geometry> geom_;
    const std::shared_ptr<Reference> ref_;
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

  public:
    MOFile(const std::shared_ptr<Geometry>, const std::shared_ptr<Reference>);
    ~MOFile();

    const std::shared_ptr<Geometry> geom() const { return geom_; };

    // creates integral files and returns the core energy.
    double create_Jiiii(const int, const int);
    const int sizeij() const { return sizeij_; };
    double mo1e(const size_t i) const { return mo1e_[i]; };
    double mo2e(const size_t i, const size_t j) const { return mo2e_[i+j*sizeij_]; };
    // strictly i <= j, k <= l
    double mo2e(const int i, const int j, const int k, const int l) const { return mo2e(address_(i,j), address_(k,l)); }; 
    double mo1e(const int i, const int j) const { return mo1e(address_(i,j)); }; 
    double* core_fock_ptr() { return core_fock_.get(); };
    double* mo1e_ptr() { return mo1e_.get(); };
    double* mo2e_ptr() { return mo2e_.get(); };

    std::shared_ptr<DF_Half> mo2e_1ext() { return mo2e_1ext_; };
    const double* const mo2e_1ext_ptr() const { return mo2e_1ext_->data(); };
    void update_1ext_ints(const std::vector<double>& coeff);

};

#endif
