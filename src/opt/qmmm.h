//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: qmmm.h
// Copyright (C) 2017 Toru Shiozaki
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


#ifndef __SRC_OPT_QMMM_H
#define __SRC_OPT_QMMM_H

#include <functional>
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/wfn/get_energy.h>
#include <src/opt/constraint.h>
#include <src/util/muffle.h>

namespace bagel {

class QMMM {
  public:
    QMMM() { }
    virtual ~QMMM() { }

    virtual void edit_input(std::shared_ptr<const Geometry> current) const = 0;
    virtual std::tuple<double,std::shared_ptr<GradFile>> do_grad(const int natom) const = 0;
};


class QMMM_Tinker : public QMMM {
  private:
    void edit_tinker_input(std::shared_ptr<const Geometry> current) const;

  public:
    QMMM_Tinker() : QMMM() { }

    void edit_input(std::shared_ptr<const Geometry> current) const override;
    std::tuple<double,std::shared_ptr<GradFile>> do_grad(const int natom) const override;
};

}
#endif
