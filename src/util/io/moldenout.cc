//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenout.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/util/io/moldenout.h>
#include <src/util/atommap.h>
#include <src/wfn/zreference.h>
#include <src/wfn/relreference.h>

using namespace bagel;
using namespace std;

/************************************************************************************
************************************************************************************/

MoldenOut::MoldenOut(string filename) : MoldenIO(filename), ofs_(filename) {
  ofs_ << "[Molden Format]" << endl;
}

MoldenOut& MoldenOut::operator<< (shared_ptr<const Molecule> mol) {
  mol_ = mol;
  write_geom();
  write_aos();

  return *this;
}

MoldenOut& MoldenOut::operator<< (shared_ptr<const Reference> ref) {
  ref_ = ref;

  auto zref = dynamic_pointer_cast<const ZReference>(ref_);
  auto rref = dynamic_pointer_cast<const RelReference>(ref_);
  if (!zref && !rref) {
    if (ref->coeff())
      write_mos();
  } else if (zref) {
    if (zref->zcoeff())
      write_mos_complex();
  } else if (rref) {
    if (rref->relcoeff_full())
      write_mos_relativistic();
  }

  if (!ref->prop_freq().empty())
    write_freq();

  return *this;
}

void MoldenOut::write_geom() {
  const int num_atoms = mol_->natom();
  stringstream ss;

  ss << "[Atoms] Angs" << endl;

  for (int i = 0; i < num_atoms; ++i) {
    shared_ptr<const Atom> cur_atom = mol_->atoms(i);

    const string cur_name = cur_atom->name();
    ss << setw(2) << cur_name << setw(8)  << i+1;

    if (to_lower(cur_name) != "q")
      ss << fixed << " " << setw(20) << cur_atom->atom_number();
    else
      ss << fixed << " " << setw(20) << setprecision(12) << cur_atom->atom_charge(); 
    ss << scientific;

    const array<double,3> cur_pos = cur_atom->position();

    ss << setw(20) << setprecision(12) << cur_pos[0]*au2angstrom__
       << setw(20) << setprecision(12) << cur_pos[1]*au2angstrom__
       << setw(20) << setprecision(12) << cur_pos[2]*au2angstrom__ << endl;
  }

  ofs_ << ss.str();
}


void MoldenOut::write_aos() {
  vector<shared_ptr<const Atom>> atoms = mol_->atoms();
  const int num_atoms = mol_->natom();
  stringstream ss;

  /************************************************************
  *  Print GTO section                                        *
  ************************************************************/
  ss << "[GTO]" << endl;

  AtomMap am;
  auto iatom = atoms.begin();
  for (int ii = 0; ii != num_atoms; ++iatom, ++ii) {
    ss << ii+1 << endl;

    vector<shared_ptr<const Shell>> shells = (*iatom)->shells();
    for (auto& ishell : shells) {
      string ang_l = am.angular_string(ishell->angular_number());
      vector<double> exponents = ishell->exponents();

      int num_contracted = ishell->contractions().size();
      for (int jj = 0; jj < num_contracted; ++jj) {
        pair<int,int> range = ishell->contraction_ranges(jj);

        ss << setw(2) << ang_l << setw(12) << range.second - range.first << endl;
        for (int kk = range.first; kk < range.second; ++kk) {
          ss << setiosflags(ios_base::scientific)
             << setw(24) << setprecision(12) << exponents[kk]
             << setw(24) << setprecision(12)
             << ishell->contractions(jj)[kk]*denormalize(ishell->angular_number(), exponents[kk]) << endl;
        }
      }
    }
    ss << endl;
  }
  ss << endl;
  ofs_ << ss.str();
}


void MoldenOut::write_mos() {
  stringstream ss;
  if (mol_->spherical())
    ss << "[5D]" << endl << "[7F]" << endl << "[9G]" << endl;
  ss << "[MO]" << endl;

  const int num_mos = !ref_->coeffB() ? ref_->coeff()->mdim() : ref_->coeffB()->mdim();
  const VectorB eig = ref_->eig();
  const VectorB occup = ref_->occup();
  const VectorB eigB = ref_->eigB();
  const VectorB occupB = ref_->occupB();

  for (int i = 0; i != num_mos; ++i) {
    // alpha
    ss << " Ene=" << setw(12) << setprecision(6) << fixed << (eig.empty() ? 0.0 : eig[i]) << endl;
    ss << " Spin=" << "  Alpha" << endl;
    ss << " Occup=" << setw(12) << occup[i] << endl;
    write_mo_one(ss, !ref_->coeffB() ? ref_->coeff()->element_ptr(0, i) : ref_->coeffA()->element_ptr(0, i));
    // possibly beta
    if (ref_->coeffB()) {
      ss << " Ene=" << setw(12) << setprecision(6) << fixed << (eigB.empty() ? 0.0 : eigB[i]) << endl;
      ss << " Spin=" << "  Beta" << endl;
      ss << " Occup=" << setw(12) << occupB[i] << endl;
      write_mo_one(ss, ref_->coeffB()->element_ptr(0, i));
    }
  }
  ofs_ << ss.str();
}


void MoldenOut::write_mos_complex() {
  stringstream ss;
  auto zref = dynamic_pointer_cast<const ZReference>(ref_);
  if (!zref)
    throw logic_error("write_mos_complex should be called with ZReference"); 

  if (mol_->spherical())
    ss << "[5D]" << endl << "[7F]" << endl << "[9G]" << endl;
  ss << "[MO]" << endl;

  const int num_mos = zref->zcoeff()->mdim();
  const VectorB eig = ref_->eig();
  const VectorB occup = ref_->occup();
  const VectorB eigB = ref_->eigB();
  const VectorB occupB = ref_->occupB();

  for (int i = 0; i != num_mos; ++i) {
    // Caution! complex coeff is only implemented for spin-free orbitals
    ss << " Ene=" << setw(12) << setprecision(6) << fixed << (eig.empty() ? 0.0 : eig[i]) << endl;
    ss << " Spin=" << "  Alpha" << endl;
    ss << " Occup=" << setw(12) << occup[i] << endl;
    write_mo_one(ss, zref->zcoeff()->element_ptr(0, i));
  }
  ofs_ << ss.str();
}


void MoldenOut::write_mos_relativistic() {
  stringstream ss;
  auto relref = dynamic_pointer_cast<const RelReference>(ref_);
  if (!relref)
    throw logic_error("write_mos_relativistic should be called with RelReference"); 

  if (mol_->spherical())
    ss << "[5D]" << endl << "[7F]" << endl << "[9G]" << endl;
  ss << "[MO]" << endl;

  const int num_mos = relref->relcoeff_full()->mdim()/2;
  const int num_aos = relref->relcoeff_full()->ndim()/4;
  const VectorB eig = ref_->eig();
  const VectorB occup = ref_->occup();
  const VectorB eigB = ref_->eigB();
  const VectorB occupB = ref_->occupB();

  auto prep = [num_aos](const complex<double>* input) {
    vector<molden_impl::complex4> out(num_aos);
    for (int i = 0; i != 4; ++i) {
      for (int j = 0; j != num_aos; ++j)
        out[j].data[i] = *input++;
    }
    return out;
  };

  for (int i = 0; i != num_mos; ++i) {
    {
      ss << " Ene=" << setw(12) << setprecision(6) << fixed << (eig.empty() ? 0.0 : eig[i]) << endl;
      ss << " Spin=" << "  Alpha" << endl;
      ss << " Occup=" << setw(12) << occup[i] << endl;
      const vector<molden_impl::complex4> tmp = prep(relref->relcoeff_full()->element_ptr(0, i*2)); 
      write_mo_one(ss, tmp.data());
    } {
      ss << " Ene=" << setw(12) << setprecision(6) << fixed << (eigB.empty() ? 0.0 : eigB[i]) << endl;
      ss << " Spin=" << "  Beta" << endl;
      ss << " Occup=" << setw(12) << occupB[i] << endl;
      const vector<molden_impl::complex4> tmp = prep(relref->relcoeff_full()->element_ptr(0, i*2+1)); 
      write_mo_one(ss, tmp.data());
    }
  }
  ofs_ << ss.str();
}


void MoldenOut::write_freq() {
  stringstream ss;
  const int num_atoms = mol_->natom();
 /************************************************************
  *  Print FREQ section                                      *
  ************************************************************/
  ss << "[FREQ]" << endl;
  for (int i = 0; i < 3*num_atoms; ++i)
    ss << setw(16) << setprecision(10) << ref_->prop_freq(i) << endl;

  ss << "[INT]" << endl;
  for (int i = 0; i < 3*num_atoms; ++i)
    ss << setw(16) << setprecision(10) << ref_->prop_ir(i) << endl;

  ss << "[FR-COORD]" << endl;
  for (int i = 0; i < num_atoms; ++i) {
    shared_ptr<const Atom> cur_atom = ref_->geom()->atoms(i);

    const string cur_name = cur_atom->name();
    const array<double,3> cur_pos = cur_atom->position();

    // molden format specifies that this section is in Bohr
    ss << setw(2) << cur_name << setw(20) << setprecision(12) << cur_pos[0]
                              << setw(20) << setprecision(12) << cur_pos[1]
                              << setw(20) << setprecision(12) << cur_pos[2] << endl;
  }

  ss << "[FR-NORM-COORD]" << endl;

  for (int i = 0; i <3*num_atoms; ++i) { // number of normal modes
    ss << setw(20) <<  "vibration    " <<  i+1 << endl;
    for (int j = 0; j < num_atoms; ++j) { // number of atoms
      ss << setw(20) << setprecision(12) << ref_->prop_eig()->element(j*3+0,i)
         << setw(20) << setprecision(12) << ref_->prop_eig()->element(j*3+1,i)
         << setw(20) << setprecision(12) << ref_->prop_eig()->element(j*3+2,i) << endl;
    }
  }
  ofs_ << ss.str();
}
