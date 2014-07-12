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


#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <src/integral/rys/eribatch.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/comprys/complexsmalleribatch.h>
#include <src/integral/comprys/complexmixederibatch.h>
#include <src/integral/libint/libint.h>
#include <src/wfn/geometry_london.h>
#include <src/wfn/geometry_connect.h>
#include <src/io/moldenin.h>
#include <src/util/constants.h>
#include <src/math/quatern.h>

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


#if 0
Geometry_London::Geometry_London(const shared_ptr<const PTree> geominfo) : Geometry(geominfo) {
/*
  // members of Molecule
  spherical_ = true;
  lmax_ = 0;

  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);

  // symmetry
  symmetry_ = to_lower(geominfo->get<string>("symmetry", "c1"));

  // cartesian or not.
  const bool cart = geominfo->get<bool>("cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }
  const bool angstrom = geominfo->get<bool>("angstrom", false);

  // static external magnetic field
  magnetic_field_ = geominfo->get_array<double,3>("magnetic_field", {{0.0, 0.0, 0.0}});
  const bool tesla = geominfo->get<bool>("tesla", false);
  if (tesla) for (int i=0; i!=3; i++) magnetic_field_[i] /= au2tesla__;

  // London basis or Gaussian w/ common origin
  const string basis_type = to_lower(geominfo->get<string>("basis_type", "giao"));
  if (basis_type == "giao" || basis_type == "london") london_ = true;
  else if (basis_type == "gaussian") london_ = false;
  else throw runtime_error("Invalid basis type entered - should be london or gaussian");

  *//* Set up atoms_ *//*
  basisfile_ = to_lower(geominfo->get<string>("basis", ""));
  if (basisfile_ == "") {
    throw runtime_error("There is no basis specification");
  } else if (basisfile_ == "molden") {
    string molden_file = geominfo->get<string>("molden_file", "");
    if (molden_file == "") throw runtime_error("No molden_in file provided");

    MoldenIn mfs(molden_file, spherical_);
    mfs.read();
    mfs >> atoms_;
  } else {

    // read the default basis file
    const shared_ptr<const PTree> bdata = PTree::read_basis(basisfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_basis");

    auto atoms = geominfo->get_child("geometry");
    for (auto& a : *atoms) {
      atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem));
    }
  }
  if (atoms_.empty()) throw runtime_error("No atoms specified at all");
  for (auto& i : atoms_)
    if (symmetry_ != "c1" && i->dummy())
      throw runtime_error("External point charges are only allowed in C1 calculations so far.");

  auxfile_ = to_lower(geominfo->get<string>("df_basis", ""));  // default value for non-DF HF.
  if (!auxfile_.empty()) {
    // read the default aux basis file
    const shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (basisfile_ == "molden") {
      for(auto& iatom : atoms_) {
        if (!iatom->dummy()) {
          aux_atoms_.push_back(make_shared<const Atom>(spherical_, iatom->name(), iatom->position(), auxfile_, make_pair(auxfile_, bdata), elem));
        } else {
          // we need a dummy atom here to be consistent in gradient computations
          aux_atoms_.push_back(iatom);
        }
      }
    } else {
      auto atoms = geominfo->get_child("geometry");
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem, true));
    }
  }

  // Misc
  aux_merged_ = false;

  common_init1();

  print_atoms();

  if (nonzero_magnetic_field()) {
    cout << "  Applied magnetic field:  (" << setprecision(4) << setw(7) << magnetic_field_[0] << ", "
                                                              << setw(7) << magnetic_field_[1] << ", "
                                                              << setw(7) << magnetic_field_[2] << ") a.u." << endl;
    const double fieldsqr = magnetic_field_[0]*magnetic_field_[0] + magnetic_field_[1]*magnetic_field_[1] + magnetic_field_[2]*magnetic_field_[2];
    cout << setprecision(0) << "  Field strength = " << au2tesla__*sqrt(fieldsqr) << " T" << endl << endl;
  }

  common_init2(true, overlap_thresh_);

  // static external electric field
  external_[0] = geominfo->get<double>("ex", 0.0);
  external_[1] = geominfo->get<double>("ey", 0.0);
  external_[2] = geominfo->get<double>("ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
*/
}


void Geometry_London::common_init2(const bool print, const double thresh, const bool nodf) {

  // compute london phase factors & add them to shells
  magnetic();

  if (london_ && nonzero_magnetic_field()) cout << "  Using London orbital basis to enforce gauge-invariance" << endl;
  if (!london_ && nonzero_magnetic_field()) cout << "  Using a common gauge origin - NOT RECOMMENDED for accurate calculations.  (Use a London orbital basis instead.)" << endl;
  if (!nonzero_magnetic_field()) cout << "  Zero magnetic field - This computation would be more efficient with a Gaussian basis set." << endl;

  // symmetry set-up
  plist_ = make_shared<Petite>(atoms_, symmetry_);
  nirrep_ = plist_->nirrep();

  if (print) {
    cout << endl;
    cout << "  Number of basis functions: " << setw(8) << nbasis() << endl;
    cout << "  Number of electrons      : " << setw(8) << nele() << endl << endl;
  }

  nuclear_repulsion_ = compute_nuclear_repulsion();

  if (!auxfile_.empty() && !nodf) {
    if (print) cout << "  Number of auxiliary basis functions: " << setw(8) << naux() << endl << endl;
    cout << "  Since a DF basis is specified, we compute 2- and 3-index integrals:" << endl;
    cout << "    o Being stored without compression. Storage requirement is "
         << setprecision(3) << static_cast<size_t>(naux_)*nbasis()*nbasis()*1.6e-8 << " GB" << endl;
    Timer timer;
    compute_integrals(thresh, nodf);
    cout << "        elapsed time:  " << setw(10) << setprecision(2) << timer.tick() << " sec." << endl << endl;
  }
}
#endif

void Geometry_London::compute_integrals(const double thresh, const bool nodf) {
    df_ = form_fit<ComplexDFDist_ints<ComplexERIBatch>>(thresh, true); // true means we construct J^-1/2
}


#if 0
Geometry_London::Geometry_London(const Geometry_London& o, shared_ptr<const PTree> geominfo, const bool discard)
  : Geometry(o, geominfo, discard) {

/*
  : schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_), london_(o.london_) {

  // members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  external_ = o.external_;
  atoms_ = o.atoms_;
  aux_atoms_ = o.aux_atoms_;
  gamma_ = o.gamma_;
  magnetic_field_ = o.magnetic_field_;

  // check all the options
  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", schwarz_thresh_);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", overlap_thresh_);
  symmetry_ = to_lower(geominfo->get<string>("symmetry", symmetry_));

  // London basis or Gaussian w/ common origin
  const string basis_type = to_lower(geominfo->get<string>("basis_type", london_ ? "giao" : "gaussian"));
  if (basis_type == "giao" || basis_type == "london") london_ = true;
  else if (basis_type == "gaussian") london_ = false;
  else throw runtime_error("Invalid basis type entered - should be london or gaussian");

  spherical_ = !geominfo->get<bool>("cartesian", !spherical_);

  // check if magnetic field has changed
  auto newfield = geominfo->get_child_optional("magnetic_field");
  if (newfield) {
    magnetic_field_ = geominfo->get_array<double,3>("magnetic_field");
    const bool tesla = geominfo->get<bool>("tesla", false);
    if (tesla) for (int i=0; i!=3; i++) magnetic_field_[i] /= au2tesla__;
  }

  // check if we need to construct shells and integrals
  auto atoms = geominfo->get_child_optional("geometry");
  const string prevbasis = basisfile_;
  basisfile_ = to_lower(geominfo->get<string>("basis", basisfile_));
  // if so, construct atoms
  if (prevbasis != basisfile_ || atoms || newfield) {
    atoms_.clear();
    const shared_ptr<const PTree> bdata = PTree::read_basis(basisfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem));
    } else {
      for (auto& a : o.atoms_)
        atoms_.push_back(make_shared<const Atom>(*a, spherical_, basisfile_, make_pair(basisfile_, bdata), elem));
    }
  }
  const string prevaux = auxfile_;
  auxfile_ = to_lower(geominfo->get<string>("df_basis", auxfile_));
  if (prevaux != auxfile_ || atoms || newfield) {
    aux_atoms_.clear();
    const shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem));
    } else {
      for (auto& a : o.atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(*a, spherical_, auxfile_, make_pair(auxfile_, bdata), elem));
    }
  }

  common_init1();
  if (prevbasis != basisfile_ || prevaux != auxfile_ || atoms || newfield) {
    // discard the previous one before we compute the new one. Note that df_'s are mutable... too bad, I know..
    if (discard)
      o.discard_df();
    common_init2(true, overlap_thresh_);
  } else {
    df_ = o.df_;
    dfs_ = o.dfs_;
    dfsl_ = o.dfsl_;
    common_init2(true, overlap_thresh_, true *//* not to calculate integrals *//*);
  }
  external_ = o.external_;

  if (nonzero_magnetic_field()) {
    cout << "  Applied magnetic field:  (" << setprecision(4) << setw(7) << magnetic_field_[0] << ", "
                                                              << setw(7) << magnetic_field_[1] << ", "
                                                              << setw(7) << magnetic_field_[2] << ") a.u." << endl;
    const double fieldsqr = magnetic_field_[0]*magnetic_field_[0] + magnetic_field_[1]*magnetic_field_[1] + magnetic_field_[2]*magnetic_field_[2];
    cout << setprecision(0) << "  Field strength = " << au2tesla__*sqrt(fieldsqr) << " T" << endl << endl;
  }
*/
}


// used in SCF initial guess.
Geometry_London::Geometry_London(const vector<shared_ptr<const Atom>> atoms, const shared_ptr<const PTree> geominfo)
 : Geometry(atoms, geominfo) {

/*
  spherical_ = true;
  lmax_ = 0;

  atoms_ = atoms;

  auto atomlist = geominfo->get_child_optional("geometry");
  const bool angstrom = geominfo->get<bool>("angstrom", false);

  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);

  // cartesian or not. Look in the atoms info to find out
  spherical_ = atoms.front()->spherical();
  // basis
  auxfile_ = geominfo->get<string>("df_basis", "");
  if (!auxfile_.empty()) {
    // read the default basis file
    const shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (atomlist) {
      for (auto& i : *atomlist)
        aux_atoms_.push_back(make_shared<const Atom>(i, spherical_, angstrom, make_pair(auxfile_, bdata), elem, true));
    } else {
      // in the molden case
      for (auto& i : atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(i->spherical(), i->name(), i->position(), auxfile_, make_pair(auxfile_, bdata), elem));
    }
  }
  // symmetry
  symmetry_ = geominfo->get<string>("symmetry", "c1");

  common_init1();

  print_atoms();

  common_init2(true, overlap_thresh_);

  // static external field
  external_[0] = geominfo->get<double>("ex", 0.0);
  external_[1] = geominfo->get<double>("ey", 0.0);
  external_[2] = geominfo->get<double>("ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
*/
}
#endif

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

