//
// BAGEL - Parallel electron correlation program.
// Filename: relshell.h 
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Matthew Kelley and Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_REL_RELSHELL_H
#define __SRC_REL_RELSHELL_H

#include <array>
#include <memory>
#include <src/util/matrix.h>
#include <src/scf/shell.h>

namespace bagel {

class RelShell : public Shell {
  protected:
    std::array<std::shared_ptr<const Matrix>,3> small_; 

    const std::shared_ptr<const Shell> aux_inc_;
    const std::shared_ptr<const Shell> aux_dec_;

    std::shared_ptr<const Matrix> overlap_compute_() const;
    std::array<std::shared_ptr<const Matrix>,3> moment_compute_(const std::shared_ptr<const Matrix> overlap) const;

  public:
    RelShell(const std::shared_ptr<const Shell> o);

    const std::shared_ptr<const Matrix> small(const int i) const { return small_[i]; }

    const std::shared_ptr<const Shell> aux_inc() const { return aux_inc_; }
    const std::shared_ptr<const Shell> aux_dec() const { return aux_dec_; }
};

}

#endif
