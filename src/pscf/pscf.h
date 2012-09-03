//
// BAGEL - Parallel electron correlation program.
// Filename: pscf.h
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


#ifndef __src_pscf_pscf_h
#define __src_pscf_pscf_h

#include <cstddef>
#include <memory>
#include <src/util/pcompfile.h>
#include <src/pscf/pgeometry.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pfock.h>
#include <src/pscf/ptildex.h>
#include <src/util/pdata.h>
#include <src/pscf/pcoeff.h>
#include <src/rysint/eribatch.h>

namespace bagel {

class PSCF {
  protected:
    const std::shared_ptr<PGeometry> geom_;
    const std::shared_ptr<POverlap> overlap_;
    const std::shared_ptr<PHcore> hcore_;
    std::shared_ptr<PTildeX> tildex_;
    std::shared_ptr<PMatrix1e> aodensity_;
    std::shared_ptr<PCoeff> coeff_;
    std::vector<double> schwarz_;

    virtual void init_schwarz();
    bool direct_;

    // created at constructor, deleted at destructor.
    double* eig_;
    void print_eig(const double*) const;
    std::shared_ptr<PCompFile<ERIBatch> > ao_eri_;

    double obtain_energy(const PMatrix1e&, const PMatrix1e&, const PMatrix1e&);

  public:
    PSCF(const std::shared_ptr<PGeometry>);
    PSCF(const std::shared_ptr<PGeometry>, std::shared_ptr<PMatrix1e>);
    ~PSCF();

    virtual void compute();
    const std::shared_ptr<PGeometry> geom() { return geom_;} ;
    const std::shared_ptr<PCoeff> coeff() { return coeff_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    // Returning a copy of eig_. Returns std::vector<double>
    const std::vector<double> eig() const {
      const size_t length = (2*geom_->K()+1) * geom_->nbasis();
      const std::vector<double> out(eig_, eig_+length);
      return out;
    }
    std::shared_ptr<PCompFile<ERIBatch> > ao_eri() { return ao_eri_; };
};

}

#endif

