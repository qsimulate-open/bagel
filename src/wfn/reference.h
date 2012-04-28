//
// Newint - Parallel electron correlation program.
// Filename: reference.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef _NEWINT_WFN_REFERENCE_H
#define _NEWINT_WFN_REFERENCE_H

#include <memory>
#include <vector>
#include <src/scf/coeff.h>
#include <src/scf/hcore.h>
#include <src/scf/geometry.h>

class Reference {

  protected:
    std::shared_ptr<Geometry> geom_; 
    std::shared_ptr<Coeff> coeff_;
    std::shared_ptr<Hcore> hcore_;
    std::vector<double> schwarz_;
    std::vector<double> eig_;

    const int nclosed_;
    const int nact_;
    const int nvirt_;

  public:
    Reference(std::shared_ptr<Geometry> g, std::shared_ptr<Coeff> c,
              std::shared_ptr<Hcore> h, const std::vector<double>& s,
              const int& nclo, const int& nact, const int& nvirt);

    ~Reference() {};

    std::shared_ptr<Geometry> geom() { return geom_; };
    std::vector<double> schwarz() { return schwarz_; };
    std::shared_ptr<Hcore> hcore() { return hcore_; };
    const std::shared_ptr<Coeff> coeff() { return coeff_; };
    const std::shared_ptr<const Coeff> coeff() const { return coeff_; };
    void set_coeff(const std::shared_ptr<Coeff> c) { coeff_ = c; };

    void set_eig(const std::vector<double>& eig) { eig_ = eig; };
    std::vector<double> eig() const { return eig_; };

    int nclosed() const { return nclosed_; };
    int nact() const { return nact_; };
    int nvirt() const { return nvirt_; };

};

#endif
