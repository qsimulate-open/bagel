//
// BAGEL - Parallel electron correlation program.
// Filename: reference_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/london/reference_london.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Reference_London)

using namespace std;
using namespace bagel;

Reference_London::Reference_London(shared_ptr<const Geometry_London> g, shared_ptr<const ZCoeff> c,
                                   const int _nclosed, const int _nact, const int _nvirt,
                                   const double en) : Reference() {
  this->cgeom_ = g;
  this->geom_ = nullptr;
  this->energy_ = en;
  this->zhcore_ = make_shared<ZHcore>(cgeom_);
  this->hcore_ = nullptr;
  this->nclosed_ = _nclosed;
  this->nact_ = _nact;
  this->nvirt_ = _nvirt;
  this->nstate_ = 1;
  this->ciwfn_ = nullptr;
  this->rdm1_ = std::vector<std::shared_ptr<RDM<1>>>();
  this->rdm2_ = std::vector<std::shared_ptr<RDM<2>>>();
  this->rdm1_av_ = nullptr;
  this->rdm2_av_ = nullptr;

  // we need to make sure that all the quantities are consistent in every MPI process
  if (c) {
    mpi__->broadcast(const_pointer_cast<ZCoeff>(c)->data(), c->size(), 0);
    zcoeff_ = c;
    coeff_ = nullptr;
  }
}


std::shared_ptr<Reference> Reference_London::project_coeff(std::shared_ptr<const Geometry> geomin) const {
  throw std::logic_error("Reference_London::project_coeff(...) should not be called");
}
