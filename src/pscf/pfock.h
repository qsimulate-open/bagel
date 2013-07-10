//
// BAGEL - Parallel electron correlation program.
// Filename: pfock.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __src_pscf_pfock_h
#define __src_pscf_pfock_h

#include <src/pscf/pgeometry.h>
#include <src/pscf/pmatrix1e.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcompfile.h>
#include <src/integral/rys/eribatch.h>
#include <memory>
#include <vector>
#include <complex>

namespace bagel {

class PFock : public PMatrix1e {
  protected:
    const std::shared_ptr<PFock> previous_;
    const std::shared_ptr<PMatrix1e> density_;
    void pfock_two_electron_part();

    std::vector<double> schwarz_;
    int S2_;
    bool direct_;

    std::shared_ptr<PCompFile<ERIBatch>> file_;

  public:
    PFock(const std::shared_ptr<PGeometry>, const std::shared_ptr<PFock>,
          const std::shared_ptr<PMatrix1e>, const std::vector<double>&, const int, const bool dir);
    PFock(const std::shared_ptr<PGeometry>, const std::shared_ptr<PFock>,
          const std::shared_ptr<PMatrix1e>, const std::vector<double>&, const int, const bool dir, std::shared_ptr<PCompFile<ERIBatch>>);
    PFock(const std::shared_ptr<PGeometry>, const std::shared_ptr<PHcore>);
    ~PFock();

    bool direct() const { return direct_; };

};

}

#endif

