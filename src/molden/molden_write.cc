//
// Newint - Parallel electron correlation program.
// Filename: molden_write.cc
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <map>

#include <src/scf/geometry.h>
#include <src/molden/molden.h>
#include <src/scf/atom.h>

using namespace std;

/************************************************************************************
*  write_geo( shared_ptr<Geometry> geo, const string molden_file )                      *
*                                                                                   *
*  Writes a molden file. Just the geometry though (Atoms section)                   *
*                                                                                   *
*  TODO: Check whether this actually works                                          *
************************************************************************************/
void Molden::write_geo(const shared_ptr<Geometry> geo, const string molden_file) {
   const int num_atoms = geo->natom();

   ofstream m_out;
   m_out.open(molden_file);

   m_out << "[Molden Format]" << endl;

   m_out << "[Atoms] Angs" << endl;

   for(int i = 0; i < num_atoms; ++i) {
      shared_ptr<Atom> cur_atom = geo->atoms(i);
      
      const string cur_name = cur_atom->name();
      const int cur_number = cur_atom->atom_number();
      const vector<double> cur_pos = cur_atom->position();

      m_out << setw(8) << cur_name << setw(8)  << i+1 
                                  << setw(8)  << cur_number 
                                  << setw(12) << setprecision(8) << cur_pos[0]/ang2bohr
                                  << setw(12) << setprecision(8) << cur_pos[1]/ang2bohr
                                  << setw(12) << setprecision(8) << cur_pos[2]/ang2bohr << endl;
   }

   m_out.close();
}
