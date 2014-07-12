//
// BAGEL - Parallel electron correlation program.
// Filename: geometry_london.cc
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


#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/comprys/complexsmalleribatch.h>
#include <src/integral/comprys/complexmixederibatch.h>
#include <src/wfn/geometry_london.h>

using namespace std;
using namespace bagel;


BOOST_CLASS_EXPORT_IMPLEMENT(Geometry_London)

Geometry_London::Geometry_London(const shared_ptr<const PTree> geominfo)
 : Geometry_base(geominfo) {

  // London basis or Gaussian w/ common origin
  const string basis_type = to_lower(geominfo->get<string>("basis_type", "giao"));
  if (basis_type == "giao" || basis_type == "london") london_ = true;
  else if (basis_type == "gaussian") london_ = false;
  else throw runtime_error("Invalid basis type entered - should be london or gaussian");

  common_init2(true, overlap_thresh_);
  get_electric_field(geominfo);
}


Geometry_London::Geometry_London(const vector<shared_ptr<const Atom>> atoms, const shared_ptr<const PTree> o)
 : Geometry_base(atoms, o) {
  common_init2(true, overlap_thresh_);
  get_electric_field(o);
}


Geometry_London::Geometry_London(const Geometry_London& o, const shared_ptr<const PTree> idata, const bool discard)
 : Geometry_base(o, idata, discard), london_(o.london_) {

  // London basis or Gaussian w/ common origin
  const string basis_type = to_lower(idata->get<string>("basis_type", london_ ? "giao" : "gaussian"));
  if (basis_type == "giao" || basis_type == "london") london_ = true;
  else if (basis_type == "gaussian") london_ = false;
  else throw runtime_error("Invalid basis type entered - should be london or gaussian");

  auto atoms = idata->get_child_optional("geometry");
  auto newfield = idata->get_child_optional("magnetic_field");
  if (o.basisfile_ != basisfile_ || o.auxfile_ != auxfile_ || atoms || newfield) {
    // discard the previous one before we compute the new one. Note that df_'s are mutable... too bad, I know..
    if (discard)
      o.discard_df();
    common_init2(true, overlap_thresh_);
  } else {
    df_ = o.df_;
    dfs_ = o.dfs_;
    dfsl_ = o.dfsl_;
    common_init2(true, overlap_thresh_, true /* not to calculate integrals */);
  }
  external_ = o.external_;
}


void Geometry_London::compute_integrals(const double thresh, const bool nodf) {
    df_ = form_fit<ComplexDFDist_ints<ComplexERIBatch>>(thresh, true); // true means we construct J^-1/2
}


shared_ptr<const Geometry_London> Geometry_London::relativistic(const bool do_gaunt) const {
  cout << "  *** Geometry_London (Relativistic) ***" << endl;
  Timer timer;
  // basically the same
  auto geom = make_shared<Geometry_London>(*this);

  // except for atoms_->shells
  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->relativistic(magnetic_field_, london_));
  geom->atoms_ = atom;

  geom->compute_relativistic_integrals(do_gaunt);

  cout << endl;
  timer.tick_print("Geometry_London relativistic (total)");
  cout << endl;
  return geom;
}


void Geometry_London::compute_relativistic_integrals(const bool do_gaunt) {
  df_->average_3index();
  dfs_  = form_fit<ComplexDFDist_ints<ComplexSmallERIBatch>>(overlap_thresh_, true, 0.0, true);
  if (do_gaunt) {
    dfsl_ = form_fit<ComplexDFDist_ints<ComplexMixedERIBatch>>(overlap_thresh_, true, 0.0, true);
  }

  // suppress some of the printing
  resources__->proc()->set_print_level(2);
}


void Geometry_London::discard_relativistic() const {
  dfs_.reset();
  dfsl_.reset();
}


void Geometry_London::magnetic() {

  Timer timer;
  const array<double,3> fieldin = london_ ? magnetic_field_ : array<double,3>{{0.0, 0.0, 0.0}};

  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->apply_magnetic_field(fieldin));
  atoms_ = atom;

  cout << endl;
  timer.tick_print("Magnetic field overhead");
  cout << endl;
}

