//
// Newint - Parallel electron correlation program.
// Filename: molden_read.cc
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

#include <src/molden/molden.h>
#include <src/scf/geometry.h>
#include <src/scf/atom.h>
#include <src/util/constants.h>

using namespace std;

/************************************************************************************
*  read_geo( const string molden_file )                                             *
*                                                                                   *
*  Reads a molden file and creates a bunch of shared_ptr's to Atoms                 *
*                                                                                   *
*  TODO: Lots of clean up. There are some unused things in here that may be used    *
*     in future functions.                                                          *
************************************************************************************/
vector<shared_ptr<Atom> > Molden::read_geo(const string molden_file) {
   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   int num_atoms = 0;

   /* Atom positions */
   vector< vector<double> > positions;
   /* Atom names */
   vector< string > names;
   /* Map atom number to basis info */
   map<int, vector<tuple<string, vector<double>, vector<double> > > > basis_info;
   /* spherical */
   bool is_spherical = false;
   /* Matrix of coefficients... not yet fully implemented */
   vector< vector<double> > mo_coefficients;

   /************************************************************
   *  Set up "global" regular expressions                      *
   ************************************************************/
   boost::regex gto_re("\\[GTO\\]");
   boost::regex atoms_re("\\[Atoms\\]");
   boost::regex mo_re("\\[MO\\]");
   boost::regex other_re("\\[\\w+\\]");

   boost::cmatch matches;

   /************************************************************
   *  Booleans to check and make sure each important section   *
   *  was found                                                *
   ************************************************************/
   bool found_atoms = false;
   bool found_gto = false;
   bool found_mo = false;

   double scale;

   /************************************************************
   *  At the moment, I'm not really sure where the [5D]        *
   *  keyword will be found, so I'm just going to look all     *
   *  the way through the file and then close it               *
   ************************************************************/

   string line; // Contains the current line of the file

   ifstream sph_input;
   sph_input.open(molden_file.c_str());
   if(!sph_input.is_open()){
      throw runtime_error("Molden input file not found");
   }   
   else {
      boost::regex _5d_re("\\[5[Dd]\\]");
      boost::regex _5d7f_re("\\[5[Dd]7[Ff]\\]");
      while (!sph_input.eof()) {
         getline(sph_input, line);
         if(boost::regex_search(line,_5d_re)){
            is_spherical = true;
         }
         else if(boost::regex_search(line,_5d7f_re)) {
            is_spherical = true;
         }
      }
   }
   sph_input.close();

   /************************************************************
   *  Open input stream                                        *
   ************************************************************/

   ifstream ifs;
   ifs.open(molden_file.c_str());
   if(!ifs.is_open()){
      throw runtime_error("Molden input file not found");
   }
   else {
      getline(ifs, line);

      while (!ifs.eof()){
         if (boost::regex_search(line,atoms_re)) {
            boost::regex ang_re("Angs");
            boost::regex atoms_line("(\\w{1,2})\\s+\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");

            scale = boost::regex_search(line, ang_re) ? ang2bohr : 1.0;

            getline(ifs, line);
            while(!boost::regex_search(line, other_re)){
               if (ifs.eof()) { break; }

               if(boost::regex_search(line.c_str(), matches, atoms_line)) {
                  ++num_atoms;

                  string nm(matches[1].first, matches[1].second);
                  names.push_back(nm);
                  
                  string x_str(matches[2].first, matches[2].second);
                  string y_str(matches[3].first, matches[3].second);
                  string z_str(matches[4].first, matches[4].second);

                  vector<double> pos;

                  pos.push_back(boost::lexical_cast<double>(x_str)*scale);
                  pos.push_back(boost::lexical_cast<double>(y_str)*scale);
                  pos.push_back(boost::lexical_cast<double>(z_str)*scale);
                  
                  positions.push_back(pos);

                  getline(ifs,line);
               }

               else { getline(ifs,line); }
            }

            found_atoms = true;
         }
         else if (boost::regex_search(line,gto_re)){
            getline(ifs, line);

            boost::regex atom_line("(\\d+)\\s+\\d?");
            boost::regex shell_line("([spd])\\s+(\\d+)\\s+\\S*");
            boost::regex exp_line("(\\S+)\\s+(\\S+)");
            boost::regex Dd("[Dd]");

            while(!boost::regex_search(line,other_re)){
               if (ifs.eof()) { break; }

               /* This line should be a new atom */
               if(!boost::regex_search(line.c_str(), matches, atom_line)) {
                  getline(ifs,line); continue;
               }
               vector<tuple<string,vector<double>,vector<double> > > atom_basis_info;

               string atom_no_str(matches[1].first, matches[1].second);
               int atom_no = boost::lexical_cast<int>(atom_no_str.c_str());

               string ang_l;

               getline(ifs, line);

               while(boost::regex_search(line.c_str(), matches, shell_line)) {
                  /* Now it should be a new angular shell */
                  string ang_l(matches[1].first, matches[1].second);

                  // TODO: figure out how to read and include [5D] keyword

                  string num_exp_string(matches[2].first, matches[2].second);
                  int num_exponents = boost::lexical_cast<int>(num_exp_string);

                  vector<double> exponents;
                  vector<double> coefficients;

                  for (int i = 0; i < num_exponents; ++i) {
                     getline(ifs,line);

                     boost::regex_search(line.c_str(), matches, exp_line);

                     string exp_string(matches[1].first, matches[1].second);
                     string coeff_string(matches[2].first, matches[2].second);

                     exp_string = boost::regex_replace(exp_string, Dd, "E");
                     coeff_string = boost::regex_replace(coeff_string, Dd, "E");

                     double exponent = boost::lexical_cast<double>(exp_string);
                     double coeff = boost::lexical_cast<double>(coeff_string);

                     exponents.push_back(exponent);
                     coefficients.push_back(coeff);
                  }

                  /************************************************************
                  *  Right now I have evec and cvec which have the newly      *
                  *  exponents and coefficients                               *
                  ************************************************************/
                  atom_basis_info.push_back(make_tuple(ang_l, exponents, coefficients));

                  getline(ifs, line); 
               }
               //atom_basis_info to basis_info
               basis_info.insert(pair<int,vector<tuple<string,vector<double>,vector<double> > > >(atom_no, atom_basis_info));
            }

            found_gto = true;
         }
         #if 0
         else if (boost::regex_search(line,mo_re)) {
            break; // this section needs work
            if(!found_gto) {
               cout << "GTO section hasn't been found yet!" << endl; break;
            }

            boost::regex ene_re("Ene=\\s+(\\S+)");
            boost::regex spin_re("Spin=\\s+(\\w+)");
            boost::regex occup_re("Occup=\\s+(\\S+)");
            boost::regex coeff_re("\\d+\\s+(\\S+)");

            while(!boost::regex_search(line,other_re)) {
               if (ifs.eof()) { break; }

               vector<double> movec;

               /* MO Energy */
               getline(ifs,line);
               if(!boost::regex_search(line.c_str(), matches, ene_re)) {
                  cout << "Whoops, there's a problem in with the MO section" << endl;
               }

               /* Spin */
               getline(ifs,line);
               if(!boost::regex_search(line.c_str(), matches, spin_re)) {
                  cout << "Whoops, there's a problem in with the MO section" << endl;
               }
               
               /* Occupation */
               getline(ifs,line);
               if(!boost::regex_search(line.c_str(), matches, occup_re)) {
                  cout << "Whoops, there's a problem in with the MO section" << endl;
               }

               /* Read in MO coefficients */
               for (int i = 0; i < num_basis; ++i) {
                  getline(ifs,line);
                  boost::regex_search(line.c_str(), matches, coeff_re);
                  string mo_string(matches[1].first, matches[1].second);
                  double coeff = boost::lexical_cast<double>(mo_string);
                  cout << i+1 << "   " << coeff << endl;
                  movec.push_back(coeff);
               }
               mo_coefficients.push_back(movec);
            }
         }
         #endif
         else {
            getline(ifs,line);
         }
      }
   }
   ifs.close();

   /************************************************************
   *  Check to make sure all the necessary information was     *
   *  found.                                                   *
   ************************************************************/
   if( !(found_atoms && found_gto) ){
      string message("Section not found in Molden file: ");
      if (!found_atoms){ message += "atoms "; }
      if (!found_gto)  { message += "GTO"; }
      throw runtime_error(message);
   }

   /************************************************************
   *  Now all the information is collected, it just has to be  *
   *  organized                                                *
   ************************************************************/
   vector<shared_ptr<Atom> > all_atoms;

   /************************************************************
   *  All the information is collected, now to organize it     *
   *  into atoms                                               *
   ************************************************************/

   /* Assuming the names and positions vectors are in the right order */
   vector<string>::iterator niter = names.begin();
   vector<vector<double> >::iterator piter = positions.begin();
   for(int i = 0; i < num_atoms; ++i, ++niter, ++piter){
      vector<tuple<string, vector<double>, vector<double> > > binfo = basis_info.find(i+1)->second;
      if (i == num_atoms) {
         throw runtime_error("It appears an atom was missing in the GTO section. Check your file");
      }

      /* For each atom, I need to make an atom object and stick it into a vector */
      shared_ptr<Atom> this_atom(new Atom(is_spherical, *niter, *piter, binfo));

      all_atoms.push_back(this_atom);
   }

   return all_atoms;
}
