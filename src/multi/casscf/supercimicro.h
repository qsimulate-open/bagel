//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: supercimicro.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_CASSCF_SUPERCIMICRO_H
#define __SRC_CASSCF_SUPERCIMICRO_H

#include <src/multi/casscf/casscf.h>

namespace bagel {

class SuperCIMicro {
  protected:
    std::shared_ptr<const CASSCF> casscf_;

    // input
    std::shared_ptr<const RotFile> grad_;
    std::shared_ptr<const RotFile> denom_;

    std::shared_ptr<const Matrix> fock_;
    std::shared_ptr<const Matrix> fockact_;
    std::shared_ptr<const Matrix> fockactp_;
    std::shared_ptr<const Matrix> gaa_;

    // results will be set to here
    std::shared_ptr<const RotFile> cc_;

    std::shared_ptr<RotFile> form_sigma(std::shared_ptr<const RotFile> cc) const;

  private:
    void sigma_at_at_(std::shared_ptr<const RotFile>, std::shared_ptr<RotFile>) const;
    void sigma_ai_ai_(std::shared_ptr<const RotFile>, std::shared_ptr<RotFile>) const;
    void sigma_at_ai_(std::shared_ptr<const RotFile>, std::shared_ptr<RotFile>) const;
    void sigma_ai_ti_(std::shared_ptr<const RotFile>, std::shared_ptr<RotFile>) const;
    void sigma_ti_ti_(std::shared_ptr<const RotFile>, std::shared_ptr<RotFile>) const;

  public:
    SuperCIMicro(std::shared_ptr<const CASSCF> cas, std::shared_ptr<const RotFile> g, std::shared_ptr<const RotFile> d,
                 std::shared_ptr<const Matrix> f1, std::shared_ptr<const Matrix> f2, std::shared_ptr<const Matrix> f3, std::shared_ptr<const Matrix> ga)
      : casscf_(cas), grad_(g), denom_(d), fock_(f1), fockact_(f2), fockactp_(f3), gaa_(ga) { }

    void compute();

    std::shared_ptr<const RotFile> cc() const { assert(cc_); return cc_; }
};

}

#endif
