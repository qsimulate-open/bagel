//
// Newint - Parallel electron correlation program.
// Filename: pscf_disk.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __src_pscf_pscf_disk_h
#define __src_pscf_pscf_disk_h

#include <memory>
#include <src/pscf/pscf.h>
#include <src/util/pmatrix1e.h>
#include <src/util/pcompfile.h>

namespace bagel {

class PSCF_DISK : public PSCF {
  protected:
    void store_ERI();

  public:
    PSCF_DISK(const std::shared_ptr<PGeometry>);
    PSCF_DISK(const std::shared_ptr<PGeometry>, std::shared_ptr<PMatrix1e>);
    ~PSCF_DISK();

};

}

#endif
