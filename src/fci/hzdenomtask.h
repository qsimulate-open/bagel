//
// BAGEL - Parallel electron correlation program.
// Filename: hzdenomtask.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#ifndef __SRC_FCI_HZ_DENOMTASK_H
#define __SRC_FCI_HZ_DENOMTASK_H

#include <bitset>
#include <memory>
#include <src/util/constants.h>
#include <src/fci/determinants.h>

namespace bagel {

class HZDenomTask {
  protected:
    double* out_;
    const std::bitset<nbit__> a_;
    const std::shared_ptr<const Determinants> det_;
    const double* jop_;
    const double* kop_;
    const double* h_;

  public:
    HZDenomTask(double* o, const std::bitset<nbit__>& a, const std::shared_ptr<const Determinants>& det, const double* j, const double* k, const double* diag)
      : out_(o), a_(a), det_(det), jop_(j), kop_(k), h_(diag) { }

    void compute() {
      const int nspin = det_->nspin();
      const int nspin2 = nspin*nspin;
      const int norb = det_->norb();
      double* iter = out_;
      const std::bitset<nbit__> ia = a_;
      for (auto& ib : det_->stringb()) {
        const int nopen = (ia^ib).count();
        const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
        *iter = 0.0;
        for (int i = 0; i != norb; ++i) {
          const int nia = ia[i];
          const int nib = ib[i];
          const int niab = nia + nib;
          const int Ni = (nia ^ nib);
          for (int j = 0; j != i; ++j) {
            const int nja = ia[j];
            const int njb = ib[j];
            const int Nj = (nja ^ njb);
            const int addj = niab * (nja + njb);
            *iter += jop_[j+norb*i] * 2.0 * addj - kop_[j+norb*i] * (F*Ni*Nj + addj);
          }
          *iter += h_[i] * niab - kop_[i+norb*i] * 0.5 * (Ni - niab*niab);
        }
        ++iter;
      }
    }

};

}

#endif
