//
// Newint - Parallel electron correlation program.
// Filename: atom.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/atom.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <tuple>
#include <src/osint/overlapbatch.h>

#define PI 3.1415926535897932

using namespace std;

typedef std::shared_ptr<Shell> RefShell;

Atom::Atom(const Atom& old, const vector<double>& displacement) 
: spherical_(old.spherical_), name_(old.name()), atom_number_(old.atom_number()), nbasis_(old.nbasis()), lmax_(old.lmax()) {
  assert(displacement.size() == 3 && old.position().size() == 3);
  vector<double>::const_iterator diter, oiter;
  const vector<double> opos = old.position();
  for (diter = displacement.begin(), oiter = opos.begin(); diter != displacement.end(); ++diter, ++oiter) {
    position_.push_back(*diter + *oiter); 
  }

  const vector<RefShell> old_shells = old.shells();
  for(vector<RefShell>::const_iterator siter = old_shells.begin(); siter != old_shells.end(); ++siter)
    shells_.push_back((*siter)->move_atom(displacement));
} 


Atom::Atom(const Atom& old, const double* displacement) 
: spherical_(old.spherical_), name_(old.name()), atom_number_(old.atom_number()), nbasis_(old.nbasis()), lmax_(old.lmax()) {
  vector<double>::const_iterator diter, oiter;
  const vector<double> opos = old.position();
  int i = 0;
  for (oiter = opos.begin(); oiter != opos.end(); ++i, ++oiter) {
    position_.push_back(displacement[i] + *oiter); 
  }

  const vector<RefShell> old_shells = old.shells();
  for(vector<RefShell>::const_iterator siter = old_shells.begin(); siter != old_shells.end(); ++siter)
    shells_.push_back((*siter)->move_atom(displacement));
}


Atom::Atom(const bool sph, const string nm, const vector<double>& p, const string basis_file)
: spherical_(sph), name_(nm), position_(p) {

  ifstream ifs;
  string bfile = basis_file;
  transform(bfile.begin(), bfile.end(), bfile.begin(),(int (*)(int))std::tolower);
  const string filename = "basis/" + bfile + ".basis";
  bool basis_found = false;
  ifs.open(filename.c_str());

  map<string, int> atommap = atommap_.atommap;
  map<string, int>::iterator tmpiter = atommap.find(nm); 
  if (tmpiter == atommap.end()) throw runtime_error("Unknown atom specified.");
  atom_number_ = tmpiter->second; 

  // this will be used in the construction of Basis_batch
  vector<tuple<string, vector<double>, vector<vector<double> > > > basis_info;

  if (!ifs.is_open()) {
    throw std::runtime_error("Basis file not found");
  } else {
    boost::regex first_line("^\\s*([spdfghijkl]+)\\s+([0-9eE\\+\\-\\.]+)\\s+");
    boost::regex other_line("^\\s*([0-9eE\\+\\-\\.-]+)\\s+");
    boost::regex coeff_line("([0-9eE\\+\\-\\.-]+)\\s*");
    string nnm = nm;
    nnm[0] = toupper(nnm[0]); 
    boost::regex atom_line("Atom:" + nnm +"\\s*$");
 
    // Reading a file to the end of the file
    string buffer;
    while (!ifs.eof()) {
      getline(ifs, buffer);
      if (buffer.empty()) continue;
      string::const_iterator start = buffer.begin();
      string::const_iterator end = buffer.end();
      boost::smatch what;
      if (regex_search(start, end, what, atom_line)) {
        basis_found = true;
        // temporary storage of info 
        string angular_number;
        vector<double> exponents;
        vector<vector<double> > coefficients;
        while (!ifs.eof()) {
          getline(ifs, buffer);
          if (buffer.empty()) continue;

          string::const_iterator start = buffer.begin();
          string::const_iterator end = buffer.end();

          if (regex_search(start, end, what, first_line)) {
            // if something is in temporary area, add to basis_info
            if (!exponents.empty()) {
              assert(!angular_number.empty());
              basis_info.push_back(make_tuple(angular_number, exponents, coefficients)); 
              exponents.clear();
              coefficients.clear();
            }
            const string ang(what[1].first, what[1].second);
            const string exp_str(what[2].first, what[2].second);
            exponents.push_back(boost::lexical_cast<double>(exp_str));
            angular_number = ang;

            start = what[0].second;
            vector<double> tmp;
            while(regex_search(start, end, what, coeff_line)) {
              const string coeff_str(what[1].first, what[1].second);
              tmp.push_back(boost::lexical_cast<double>(coeff_str));
              start = what[0].second;
            }
            coefficients.push_back(tmp);
          } else if (regex_search(start, end, what, other_line)) {
            const string exp_str(what[1].first, what[1].second);
            exponents.push_back(boost::lexical_cast<double>(exp_str));
            
            start = what[0].second;
            vector<double> tmp;
            while(regex_search(start, end, what, coeff_line)) {
              const string coeff_str(what[1].first, what[1].second);
              tmp.push_back(boost::lexical_cast<double>(coeff_str));
              start = what[0].second;
            }
            coefficients.push_back(tmp);
          } else {
            basis_info.push_back(make_tuple(angular_number, exponents, coefficients)); 
            break;
          }
        }
        break;
      }
    }
    if (!basis_found) throw runtime_error("Basis was not found.");

    construct_shells(basis_info);
  }

  ifs.close();

  // counting the number of basis functions belonging to this atom
  nbasis_ = 0;
  for (vector<RefShell>::iterator siter = shells_.begin(); siter != shells_.end(); ++siter) {
    const int ang = (*siter)->angular_number();
    if (spherical_) {
      nbasis_ += (*siter)->num_contracted() * (2 * ang + 1); 
    } else {
      nbasis_ += (*siter)->num_contracted() * (ang + 1) * (ang + 2) / 2; 
    }
  }

//print_basis();
}

Atom::~Atom() {

}


/* NOTE : Actually I realized that the code works with the following basis format as well

Atom:Li
s        1469.0000000              0.0007660       -0.0001200 
          220.5000000              0.0058920       -0.0009230
           50.2600000              0.0296710       -0.0046890
           14.2400000              0.1091800       -0.0176820 
            4.5810000              0.2827890       -0.0489020 
            1.5800000              0.4531230       -0.0960090 
            0.5640000              0.2747740       -0.1363800 
            0.0734500              0.0097510        0.5751020
s           0.0280500              1.0000000    
p           1.5340000              0.0227840    
            0.2749000              0.1391070    
            0.0736200              0.5003750    
p           0.0240300              1.0000000    
d           0.1239000              1.0000000    

which was the reason why the third argument was a vector of a vector.

*/


// convert basis_info to vector<Shell> 
void Atom::construct_shells(vector<tuple<string, vector<double>, vector<vector<double> > > > in) {

  for (int i = 0; i != 10; ++i) { 
    vector<vector<double> > contractions;
    vector<pair<int, int> > contranges;
    vector<double> exponents;

    // previous set of exponents... 
    vector<double> previous_exp;
    
    int offset = 0;
    for (auto biter = in.begin(); biter != in.end(); ++biter) {

      // check the angular number
      if (atommap_.angular_number(get<0>(*biter)) != i) continue;
      
      // contraction coefficient matrix
      const vector<vector<double> > conts = get<2>(*biter);

      // loop over contraction coefficients
      for (int j = 0; j != conts.front().size(); ++j) {

        // picking the current contraction coefficients (for the case there are multiple coefficients specified).
        vector<double> current;
        for (auto citer = conts.begin(); citer != conts.end(); ++citer) 
          current.push_back((*citer)[j]);

        // counting the number of zeros above and below in the segmented contractions
        int zerostart = 0;
        for (auto iter = current.begin(); iter != current.end(); ++iter) {
          if (*iter == 0.0) ++zerostart;
          else break;
        }
        int zeroend = 0;
        for (auto iter = current.rbegin(); iter != current.rend(); ++iter) {
          if (*iter == 0.0) ++zeroend;
          else break;
        }

        // first make a vector with zero
        vector<double> cont2(offset, 0.0); 
        // and add the coefficients
        cont2.insert(cont2.end(), current.begin(), current.end());
        contractions.push_back(cont2);
        contranges.push_back(make_pair(offset + zerostart, offset + current.size() - zeroend));
        assert(offset + zerostart <= offset + current.size() - zeroend);
      }
      const vector<double> exp = get<1>(*biter);
// TODO we only increment when exponents are differnent from the previous one 
      if (true || exp != previous_exp) {
        exponents.insert(exponents.end(), exp.begin(), exp.end()); 
        offset += exp.size(); 
        previous_exp = exp;
      }
    }

    // this is to do with normalization
    if (!exponents.empty()) { 
      auto citer = contranges.begin();
      for (auto iter = contractions.begin(); iter != contractions.end(); ++iter, ++citer) {
        auto eiter = exponents.begin();
        double denom = 1.0;
        for (int ii = 2; ii <= i; ++ii) denom *= 2 * ii - 1;
        for (auto diter = iter->begin(); diter != iter->end(); ++diter, ++eiter)
          *diter *= ::pow(2.0 * *eiter / PI, 0.75) * ::pow(::sqrt(4.0 * *eiter), static_cast<double>(i)) / ::sqrt(denom);

        vector<vector<double> > cont(1, *iter);
        vector<pair<int, int> > cran(1, *citer);
        RefShell current(new Shell(spherical_, position_, i, exponents, cont, cran));
        vector<RefShell> cinp(2, current); 
        OverlapBatch coverlap(cinp);
        coverlap.compute();
        const double scal = 1.0 / ::sqrt((coverlap.data())[0]);
        for (auto diter = iter->begin(); diter != iter->end(); ++diter) *diter *= scal;
      } 

      RefShell currentbatch(new Shell(spherical_, position_, i, exponents, contractions, contranges));
      shells_.push_back(currentbatch);
      lmax_ = i;
    }

  } // end of batch loop 

}


void Atom::print_basis() const {
  for(auto iter = shells_.begin(); iter != shells_.end(); ++iter)
     cout << (*iter)->show() << endl; 
}


void Atom::print() const {
  string tmp = name_;
  tmp[0] = ::toupper(tmp[0]);
  cout << "  atom = (" << setw(2) << tmp << "," << fixed << setprecision(6) <<
      setw(14) << position_[0] << "," <<
      setw(14) << position_[1] << "," << 
      setw(14) << position_[2] <<  ");" << endl;
}

