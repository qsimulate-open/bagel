//
// BAGEL - Parallel electron correlation program.
// Filename: geometry.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
#include <src/integral/rys/smalleribatch.h>
#include <src/integral/rys/mixederibatch.h>
#include <src/integral/libint/libint.h>
#include <src/wfn/geometry.h>
#include <src/wfn/geometry_connect.h>
#include <src/io/moldenin.h>
#include <src/util/constants.h>
#include <src/math/quatern.h>

using namespace std;
using namespace bagel;


BOOST_CLASS_EXPORT_IMPLEMENT(Geometry)


Geometry::Geometry(const shared_ptr<const PTree> geominfo)
 : Geometry_base(geominfo) {
  common_init2(true, overlap_thresh_);
  get_electric_field(geominfo);
}


Geometry::Geometry(const vector<shared_ptr<const Atom>> atoms, const shared_ptr<const PTree> o)
 : Geometry_base(atoms, o) {
  common_init2(true, overlap_thresh_);
  get_electric_field(o);
}


Geometry::Geometry(const Geometry& o, const shared_ptr<const PTree> idata, const bool discard)
 : Geometry_base(o, idata, discard) {
  auto atoms = idata->get_child_optional("geometry");
  if (o.basisfile_ != basisfile_ || o.auxfile_ != auxfile_ || atoms) {
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


Geometry::Geometry(const Geometry& o, const shared_ptr<const Matrix> disp, const shared_ptr<const PTree> geominfo, const bool rotate, const bool nodf)
 : Geometry_base(o, disp, geominfo, rotate, nodf) {
  common_init2(false, overlap_thresh_, nodf);
}


Geometry::Geometry(const Geometry& o, const array<double,3> disp)
 : Geometry_base(o, disp) {
  common_init2(false, overlap_thresh_);
}


/************************************************************
*  Merge info from multiple geometries to make one          *
*  supergeometry                                            *
************************************************************/
Geometry::Geometry(vector<shared_ptr<const Geometry>> nmer) {

   // members of base classes
   spherical_ = nmer.front()->spherical_;
   symmetry_ = nmer.front()->symmetry_;
   external_ = nmer.front()->external_;
   schwarz_thresh_ = nmer.front()->schwarz_thresh_;
   overlap_thresh_ = nmer.front()->overlap_thresh_;

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

void Geometry::compute_integrals(const double thresh, const bool nodf) {
#ifdef LIBINT_INTERFACE
    df_ = form_fit<DFDist_ints<Libint>>(thresh, true); // true means we construct J^-1/2
#else
    df_ = form_fit<DFDist_ints<ERIBatch>>(thresh, true); // true means we construct J^-1/2
#endif
}


shared_ptr<const Geometry> Geometry::relativistic(const bool do_gaunt) const {
  cout << "  *** Geometry (Relativistic) ***" << endl;
  Timer timer;
  // basically the same
  auto geom = make_shared<Geometry>(*this);

  // except for atoms_->shells
  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->relativistic());
  geom->atoms_ = atom;

  geom->compute_relativistic_integrals(do_gaunt);

  cout << endl;
  timer.tick_print("Geometry relativistic (total)");
  cout << endl;
  return geom;
}


void Geometry::compute_relativistic_integrals(const bool do_gaunt) {
  df_->average_3index();
  dfs_  = form_fit<DFDist_ints<SmallERIBatch>>(overlap_thresh_, true, 0.0, true);
  if (do_gaunt)
    dfsl_ = form_fit<DFDist_ints<MixedERIBatch>>(overlap_thresh_, true, 0.0, true);

  // suppress some of the printing
  resources__->proc()->set_print_level(2);
}


void Geometry::discard_relativistic() const {
  dfs_.reset();
  dfsl_.reset();
}
