//
// BAGEL - Parallel electron correlation program.
// Filename: moldenout.cc
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/io/moldenout.h>
#include <src/util/atommap.h>

using namespace bagel;
using namespace std;

/************************************************************************************
************************************************************************************/

MoldenOut::MoldenOut(string filename) : MoldenIO(filename) {
  ofs_.open(filename);
  ofs_ << "[Molden Format]" << endl;
}

MoldenOut& MoldenOut::operator<< (shared_ptr<const Geometry> geom) {
  geom_ = geom;
  write_geom();

  return *this;
}

MoldenOut& MoldenOut::operator<< (shared_ptr<const Reference> ref) {
  ref_ = ref;

  write_mos();

  return *this;
}

void MoldenOut::write_geom() {
  const int num_atoms = geom_->natom();

  ofs_ << "[Atoms] Angs" << endl;

  for (int i = 0; i < num_atoms; ++i) {
     shared_ptr<const Atom> cur_atom = geom_->atoms(i);

     const string cur_name = cur_atom->name();
     const int cur_number = cur_atom->atom_number();
     const array<double,3> cur_pos = cur_atom->position();

     ofs_ << setw(2) << cur_name << setw(8)  << i+1
                                 << setw(8)  << cur_number << setiosflags(ios_base::scientific)
                                 << setw(20) << setprecision(12) << cur_pos[0]/ang2bohr__
                                 << setw(20) << setprecision(12) << cur_pos[1]/ang2bohr__
                                 << setw(20) << setprecision(12) << cur_pos[2]/ang2bohr__ << endl;
  }
}

void MoldenOut::write_mos() {
  vector<shared_ptr<const Atom>> atoms = geom_->atoms();
  const bool is_spherical = geom_->spherical();

  const int num_atoms = geom_->natom();

  /************************************************************
  *  Print GTO section                                        *
  ************************************************************/
  ofs_ << "[GTO]" << endl;

  AtomMap am;
  auto iatom = atoms.begin();
  for (int ii = 0; ii != num_atoms; ++iatom, ++ii) {
    ofs_ << ii+1 << endl;

    vector<shared_ptr<const Shell>> shells = (*iatom)->shells();
    for (auto& ishell : shells) {
      string ang_l = am.angular_string(ishell->angular_number());
      vector<double> exponents = ishell->exponents();

      int num_contracted = ishell->contractions().size();
      for (int jj = 0; jj < num_contracted; ++jj) {
        pair<int,int> range = ishell->contraction_ranges(jj);

        ofs_ << setw(2) << ang_l << setw(8) << range.second - range.first << endl;
        for (int kk = range.first; kk < range.second; ++kk) {
          ofs_ << setiosflags(ios_base::scientific)
               << setw(20) << setprecision(8) << exponents[kk]
               << setw(20) << setprecision(8)
               << ishell->contractions(jj)[kk]*denormalize(ishell->angular_number(), exponents[kk]) << endl;
        }
      }
    }
    ofs_ << endl;
  }
  ofs_ << endl;
  if (is_spherical) ofs_ << "[5D]" << endl;
  ofs_ << "[MO]" << endl;

  const int nbasis = ref_->coeff()->ndim();
  const int num_mos = ref_->coeff()->mdim();
  int nocc = ref_->nclosed();

  const double* modata = ref_->coeff()->data();

  vector<double> eigvec = ref_->eig();
  if (eigvec.empty()) eigvec = vector<double>(num_mos,0.0);

  int i = 0;
  auto ieig = eigvec.begin();

  for (int i = 0; i < num_mos; ++ieig, ++i) {

    ofs_ << " Ene=" << setw(12) << setprecision(6) << fixed << *ieig << endl;

    /* At the moment only thinking about RHF, so assume spin is Alpha */
    ofs_ << " Spin=" << "  Alpha" << endl;

    /* At the moment, assuming occupation can be 2 or 0. Should be fine for RHF */
    string occ_string = nocc-- > 0 ? "  2.000" : "  0.000";
    ofs_ << " Occup=" << occ_string << endl;

    int j = 0;
    for (auto& iatom : atoms) {
      for (auto& ishell : iatom->shells()) {
        for (int icont = 0; icont != ishell->num_contracted(); ++icont) {
          vector<int> corder = (is_spherical ? b2m_sph_.at(ishell->angular_number()) : b2m_cart_.at(ishell->angular_number()));
          vector<double> scales = (is_spherical ? vector<double>(corder.size(), 1.0) : scaling_.at(ishell->angular_number())) ;
          for (auto& iorder : corder) {
             ofs_ << fixed << setw(4) << ++j << setw(22) << setprecision(16) << modata[iorder] / scales.at(iorder) << endl;
          }
          modata += corder.size();
        }
      }
    }
  }
}
