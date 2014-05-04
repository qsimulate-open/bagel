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
#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/rys/smalleribatch.h>
#include <src/integral/rys/mixederibatch.h>
#include <src/integral/libint/libint.h>
#include <src/wfn/geometry_london.h>
#include <src/wfn/geometry_connect.h>
#include <src/io/moldenin.h>
#include <src/util/constants.h>
#include <src/math/quatern.h>

using namespace std;
using namespace bagel;


BOOST_CLASS_EXPORT_IMPLEMENT(Geometry_London)

Geometry_London::Geometry_London(const shared_ptr<const PTree> geominfo) {

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

  /* Set up atoms_ */
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
      atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem, magnetic_field_));
    }
  }
  if (atoms_.empty()) throw runtime_error("No atoms specified at all");
  for (auto& i : atoms_)
    if (symmetry_ != "c1" && i->dummy())
      throw runtime_error("External point charges are only allowed in C1 calculations so far.");

  /* Set up aux_atoms_ */

  const std::array<double,3> magnetic_field_aux = {{0.0, 0.0, 0.0}};
/*
#if 1
  magnetic_field_aux[0] = magnetic_field_[1];
  magnetic_field_aux[1] = magnetic_field_[2];
  magnetic_field_aux[2] = magnetic_field_[0];
#else
  magnetic_field_aux[0] = 0.0;
  magnetic_field_aux[1] = 0.0;
  magnetic_field_aux[2] = 0.0;
#endif
*/

  auxfile_ = to_lower(geominfo->get<string>("df_basis", ""));  // default value for non-DF HF.
  if (!auxfile_.empty()) {
    // read the default aux basis file
    const shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (basisfile_ == "molden") {
      for(auto& iatom : atoms_) {
        if (!iatom->dummy()) {
          aux_atoms_.push_back(make_shared<const Atom>(spherical_, iatom->name(), iatom->position(), auxfile_, make_pair(auxfile_, bdata), elem, magnetic_field_aux));
        } else {
          // we need a dummy atom here to be consistent in gradient computations
          aux_atoms_.push_back(iatom);
        }
      }
    } else {
      auto atoms = geominfo->get_child("geometry");
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem, magnetic_field_aux, true));
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
}


void Geometry_London::common_init2(const bool print, const double thresh, const bool nodf) {
  if (!nonzero_magnetic_field()) cout << "  Zero magnetic field - This computation should be more efficient with a Gaussian basis set." << endl << endl;

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


void Geometry_London::compute_integrals(const double thresh, const bool nodf) {
    df_ = form_fit<DFDist_London_ints<ComplexERIBatch>>(thresh, true); // true means we construct J^-1/2
}


// suitable for geometry updates in optimization
Geometry_London::Geometry_London(const Geometry_London& o, const shared_ptr<const Matrix> displ, const shared_ptr<const PTree> geominfo, const bool rotate, const bool nodf)
  : schwarz_thresh_(o.schwarz_thresh_), gamma_(o.gamma_) {
  throw logic_error("This constructor of Geometry_London, designed for geometry optimization, has not been implemented or verified"); }

/*
  // Members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  external_ = o.external_;

  // first construct atoms using displacements
  int iat = 0;
  for (auto i = o.atoms_.begin(), j = o.aux_atoms_.begin(); i != o.atoms_.end(); ++i, ++j, ++iat) {
    array<double,3> cdispl = {{displ->element(0,iat), displ->element(1,iat), displ->element(2,iat)}};
    atoms_.push_back(make_shared<Atom>(**i, cdispl));
    aux_atoms_.push_back(make_shared<Atom>(**j, cdispl));
  }

  // second find the unique frame.
  // (i) center of chages
  if (rotate) {
    Quatern<double> oc = o.charge_center();
    Quatern<double> mc = charge_center();
    // (2) direction of the first atom
    Quatern<double> oa = o.atoms().front()->position();
    Quatern<double> ma =   atoms().front()->position();
    Quatern<double> od = oa - oc;
    Quatern<double> md = ma - mc;
    // Quaternion that maps md to od.
    od.normalize();
    md.normalize();
    Quatern<double> op = md * od;
    op[0] = 1.0 - op[0];
    op.normalize();
    Quatern<double> opd = op.dagger();

    // first subtract mc, rotate, and then add oc
    vector<shared_ptr<const Atom>> newatoms;
    vector<shared_ptr<const Atom>> newauxatoms;
    for (auto i = atoms_.begin(), j = aux_atoms_.begin(); i != atoms_.end(); ++i, ++j) {
      assert((*i)->position() == (*j)->position());
      Quatern<double> source = (*i)->position();
      Quatern<double> target = op * (source - mc) * opd + oc;
      array<double,3> cdispl = (target - source).ijk();

      newatoms.push_back(make_shared<Atom>(**i, cdispl));
      newauxatoms.push_back(make_shared<Atom>(**j, cdispl));
    }
    atoms_ = newatoms;
    aux_atoms_ = newauxatoms;

    // (3) plane of center of charges, first and second atoms.
    if (natom() > 2) {
      assert(natom() == o.natom());
      Quatern<double> oa0 = o.atoms(0)->position();
      Quatern<double> ma0 =   atoms(0)->position();
      Quatern<double> oa1 = o.atoms(1)->position();
      Quatern<double> ma1 =   atoms(1)->position();
      mc = charge_center();
      od = (oa0 - oc) * (oa1 - oc);
      md = (ma0 - mc) * (ma1 - mc);
      od[0] = 0.0;
      md[0] = 0.0;
      od.normalize();
      md.normalize();
      op = md * od;
      op[0] = 1.0 - op[0];
      op.normalize();
      opd = op.dagger();

      newatoms.clear();
      newauxatoms.clear();
      for (auto i = atoms_.begin(), j = aux_atoms_.begin(); i != atoms_.end(); ++i, ++j) {
        assert((*i)->position() == (*j)->position());
        Quatern<double> source = (*i)->position();
        Quatern<double> target = op * (source - mc) * opd + oc;
        array<double,3> cdispl = (target - source).ijk();

        newatoms.push_back(make_shared<Atom>(**i, cdispl));
        newauxatoms.push_back(make_shared<Atom>(**j, cdispl));
      }
      atoms_ = newatoms;
      aux_atoms_ = newauxatoms;
    }
  }

  common_init1();
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);
  common_init2(false, overlap_thresh_, nodf);
}
*/


Geometry_London::Geometry_London(const Geometry_London& o, const array<double,3> displ)
  : schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_), gamma_(o.gamma_) {
  throw logic_error("This constructor of Geometry_London has not been implemented or verified"); }

/*
  // members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  external_ = o.external_;

  // first construct atoms using displacements
  for (auto& i : o.atoms_) {
    atoms_.push_back(make_shared<Atom>(*i, displ));
  }

  for (auto& j : o.aux_atoms_) {
    aux_atoms_.push_back(make_shared<Atom>(*j, displ));
  }

  common_init1();
  common_init2(false, overlap_thresh_);
}
*/


Geometry_London::Geometry_London(const Geometry_London& o, shared_ptr<const PTree> geominfo, const bool discard)
  : schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_), gamma_(o.gamma_) {

  // members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  external_ = o.external_;
  atoms_ = o.atoms_;
  aux_atoms_ = o.aux_atoms_;
  magnetic_field_ = o.magnetic_field_;

  // check all the options
  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", schwarz_thresh_);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", overlap_thresh_);
  symmetry_ = to_lower(geominfo->get<string>("symmetry", symmetry_));

  spherical_ = !geominfo->get<bool>("cartesian", !spherical_);

  // check if we need to construct shells and integrals
  auto atoms = geominfo->get_child_optional("geometry");
  const string prevbasis = basisfile_;
  basisfile_ = to_lower(geominfo->get<string>("basis", basisfile_));
  // if so, construct atoms
  if (prevbasis != basisfile_ || atoms) {
    atoms_.clear();
    const shared_ptr<const PTree> bdata = PTree::read_basis(basisfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem, magnetic_field_));
    } else {
      for (auto& a : o.atoms_)
        atoms_.push_back(make_shared<const Atom>(*a, spherical_, basisfile_, make_pair(basisfile_, bdata), elem));
      throw runtime_error("Need to set it up to initialize magnetic field in this other place, I guess.");
    }
  }
  const string prevaux = auxfile_;
  auxfile_ = to_lower(geominfo->get<string>("df_basis", auxfile_));
  if (prevaux != auxfile_ || atoms) {
    aux_atoms_.clear();
    const shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    const shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem, magnetic_field_));
    } else {
      for (auto& a : o.atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(*a, spherical_, auxfile_, make_pair(auxfile_, bdata), elem));
    }
  }

  common_init1();
  if (prevbasis != basisfile_ || prevaux != auxfile_ || atoms) {
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


/************************************************************
*  Merge info from multiple geometries to make one          *
*  supergeometry                                            *
************************************************************/
Geometry_London::Geometry_London(vector<shared_ptr<const Geometry_London>> nmer) :
   schwarz_thresh_(nmer.front()->schwarz_thresh_), overlap_thresh_(nmer.front()->overlap_thresh_)
{ throw logic_error("This constructor of Geometry_London, designed for merging several into one, has not been implemented or verified"); }
#if 0
{
   // A member of Molecule
   spherical_ = nmer.front()->spherical_;
   symmetry_ = nmer.front()->symmetry_;
   external_ = nmer.front()->external_;

   /************************************************************
   * Going down the list of protected variables, merge the     *
   * data, pick the best ones, or make sure they all match     *
   ************************************************************/
   /* spherical_ should match across the vector*/
   for(auto& inmer : nmer) {
      if (spherical_ != (inmer)->spherical_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of cartesian and spherical bases");
      }
   }

   /* symmetry_ should match across the vector*/
   for(auto& inmer : nmer) {
      if (symmetry_ != (inmer)->symmetry_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of different symmetries");
      }
   }

   /* external field would hopefully match, but for now, if it doesn't, just disable */
   for(auto& inmer : nmer) {
      if(!(equal(external_.begin(), external_.end(), inmer->external_.begin()))){
         fill(external_.begin(), external_.end(), 0.0); break;
      }
   }

   /* atoms_ and aux_atoms_ can be merged */
   vector<shared_ptr<const Atom>> new_atoms;
   vector<shared_ptr<const Atom>> new_aux_atoms;
   for(auto& inmer : nmer) {
      auto iatoms = inmer->atoms();
      auto iaux = inmer->aux_atoms();

      new_atoms.insert(new_atoms.end(), iatoms.begin(), iatoms.end());
      new_aux_atoms.insert(new_aux_atoms.end(), iaux.begin(), iaux.end());
   }
   atoms_ = new_atoms;
   aux_atoms_ = new_aux_atoms;

   basisfile_ = nmer.front()->basisfile_;
   auxfile_ = nmer.front()->auxfile();
   aux_merged_ = nmer.front()->aux_merged_;
   gamma_ = nmer.front()->gamma();

   /* Use the strictest thresholds */
   for(auto& inmer : nmer) {
      schwarz_thresh_ = min(schwarz_thresh_, inmer->schwarz_thresh_);
      overlap_thresh_ = min(overlap_thresh_, inmer->overlap_thresh_);
   }

   /* Data is merged (crossed fingers), now finish */
   common_init1();
   print_atoms();
   common_init2(true,overlap_thresh_);

   // static external field
   if (external())
   cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                          << setw(7) << external_[1] << ", "
                                                                          << setw(7) << external_[2] << ") a.u." << endl << endl;
}
#endif


// used in SCF initial guess.
Geometry_London::Geometry_London(const vector<shared_ptr<const Atom>> atoms, const shared_ptr<const PTree> geominfo) {

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
    const std::array<double,3> magnetic_field_aux = {{0.0, 0.0, 0.0}};
    if (atomlist) {
      for (auto& i : *atomlist)
        aux_atoms_.push_back(make_shared<const Atom>(i, spherical_, angstrom, make_pair(auxfile_, bdata), elem, magnetic_field_aux, true));
    } else {
      // in the molden case
      for (auto& i : atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(i->spherical(), i->name(), i->position(), auxfile_, make_pair(auxfile_, bdata), elem, magnetic_field_aux));
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
}


const shared_ptr<const Matrix> Geometry_London::compute_grad_vnuc() const {
  // the derivative of Vnuc
  auto grad = make_shared<Matrix>(3, natom());
  int i = 0;
  for (auto& a : atoms_) {
    const double ac = a->atom_charge();
    for (auto& b : atoms_) {
      if (a == b) continue;
      const array<double,3> displ = a->displ(b);
      const double c = b->atom_charge() * ac;
      const double dist = a->distance(b);
      const double dist3 = dist*dist*dist;
      grad->element(0,i) += c*displ[0]/dist3;
      grad->element(1,i) += c*displ[1]/dist3;
      grad->element(2,i) += c*displ[2]/dist3;
    }
    ++i;
  }
  return grad;
}


void Geometry_London::merge_obs_aux() {
  aux_merged_ = true;
  atoms_.insert(atoms_.end(), aux_atoms_.begin(), aux_atoms_.end());
  for (auto iter = aux_offsets_.begin(); iter != aux_offsets_.end(); ++iter) {
    for (auto citer = iter->begin(); citer != iter->end(); ++citer) {
      *citer += nbasis_;
    }
  }
  offsets_.insert(offsets_.end(), aux_offsets_.begin(), aux_offsets_.end());
  nbasis_ += naux_;
}


vector<double> Geometry_London::schwarz() const {
  vector<shared_ptr<const Shell>> basis;
  for (auto aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter) {
    const vector<shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
  }
  const int size = basis.size();

  vector<double> schwarz(size * size);
  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<const Shell> b0 = basis[i0];
    for (int i1 = i0; i1 != size; ++i1) {
      const shared_ptr<const Shell> b1 = basis[i1];

      array<shared_ptr<const Shell>,4> input = {{b1, b0, b1, b0}};
#ifdef LIBINT_INTERFACE
      Libint eribatch(input);
#else
      ERIBatch eribatch(input, 1.0);
#endif
      eribatch.compute();
      const double* eridata = eribatch.data();
      const int datasize = eribatch.data_size();
      double cmax = 0.0;
      for (int xi = 0; xi != datasize; ++xi, ++eridata) {
        const double absed = fabs(*eridata);
        if (absed > cmax) cmax = absed;
      }
      schwarz[i0 * size + i1] = cmax;
      schwarz[i1 * size + i0] = cmax;
    }
  }
  return schwarz;
}


shared_ptr<const Geometry_London> Geometry_London::relativistic(const bool do_gaunt) const {
  cout << "  *** Geometry_London (Relativistic) ***" << endl;
  Timer timer;
  // basically the same
  auto geom = make_shared<Geometry_London>(*this);

  // except for atoms_->shells
  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->relativistic());
  geom->atoms_ = atom;

  geom->compute_relativistic_integrals(do_gaunt);

  cout << endl;
  timer.tick_print("Geometry_London relativistic (total)");
  cout << endl;
  return geom;
}

// TODO Need to make ComplexSmallERIBatch and ComplexMixedERIBatch?
void Geometry_London::compute_relativistic_integrals(const bool do_gaunt) {
  df_->average_3index();
//  dfs_  = form_fit<DFDist_London_ints<SmallERIBatch>>(overlap_thresh_, true, 0.0, true);
//  if (do_gaunt)
//    dfsl_ = form_fit<DFDist_London_ints<MixedERIBatch>>(overlap_thresh_, true, 0.0, true);

  // suppress some of the printing
  resources__->proc()->set_print_level(2);
}


void Geometry_London::discard_relativistic() const {
  dfs_.reset();
  dfsl_.reset();
}
