//
// Author: Toru Shiozaki
// Date  : May 2009
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
using namespace boost;

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
    regex first_line("^\\s*([spdfghijkl]+)\\s+([0-9eE\\+\\-\\.]+)\\s+");
    regex other_line("^\\s*([0-9eE\\+\\-\\.-]+)\\s+");
    regex coeff_line("([0-9eE\\+\\-\\.-]+)\\s*");
    regex atom_line("Atom:" + nm );
 
    string buffer;
    while (!ifs.eof()) {
      getline(ifs, buffer);
      if (buffer.empty()) continue;
      string::const_iterator start = buffer.begin();
      string::const_iterator end = buffer.end();
      smatch what;
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
            exponents.push_back(lexical_cast<double>(exp_str));
            angular_number = ang;

            start = what[0].second;
            vector<double> tmp;
            while(regex_search(start, end, what, coeff_line)) {
              const string coeff_str(what[1].first, what[1].second);
              tmp.push_back(lexical_cast<double>(coeff_str));
              start = what[0].second;
            }
            coefficients.push_back(tmp);
          } else if (regex_search(start, end, what, other_line)) {
            const string exp_str(what[1].first, what[1].second);
            exponents.push_back(lexical_cast<double>(exp_str));
            
            start = what[0].second;
            vector<double> tmp;
            while(regex_search(start, end, what, coeff_line)) {
              const string coeff_str(what[1].first, what[1].second);
              tmp.push_back(lexical_cast<double>(coeff_str));
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

    // convert basis_info to vector<Shell> 
    vector<tuple<string, vector<double>, vector<vector<double> > > >::const_iterator biter; 
    for (int i = 0; i != 10; ++i) { 
      vector<vector<double> > contractions;
      vector<pair<int, int> > contranges;
      vector<double> exponents;
      
      int offset = 0;
      for (biter = basis_info.begin(); biter != basis_info.end(); ++biter) {
        map<string, int> angmap = atommap_.angmap;
        map<string, int>::iterator miter = angmap.find(get<0>(*biter));  
        if (miter == angmap.end()) throw runtime_error("Unknown angular number in a basis set file.");
        const int angular = miter->second; 
        if (angular != i) continue;
        
        const vector<vector<double> > conts = get<2>(*biter);
        for (int j = 0; j != conts.front().size(); ++j) {
          vector<double> current;
          for (vector<vector<double> >::const_iterator citer = conts.begin(); citer != conts.end(); ++citer) 
            current.push_back((*citer)[j]);
          // checking zero's above and below in the segmented contractions
          int zerostart = 0;
          for (vector<double>::const_iterator iter = current.begin(); iter != current.end(); ++iter) {
            if (*iter == 0.0) ++zerostart;
            else break;
          }
          int zeroend = 0;
          for (vector<double>::const_reverse_iterator iter = current.rbegin(); iter != current.rend(); ++iter) {
            if (*iter == 0.0) ++zeroend;
            else break;
          }
          vector<double> cont2(offset, 0.0); 
          cont2.insert(cont2.end(), current.begin(), current.end());
          contractions.push_back(cont2);
          contranges.push_back(make_pair(offset + zerostart, offset + current.size() - zeroend));
          assert(offset + zerostart <= offset + current.size() - zeroend);
        }
        const vector<double> exp = get<1>(*biter);
        exponents.insert(exponents.end(), exp.begin(), exp.end()); 
        offset += exp.size(); 
      }
      if (!exponents.empty()) { 
//  normalizing
        vector<pair<int, int> >::iterator citer = contranges.begin();
        for (vector<vector<double> >::iterator iter = contractions.begin(); iter != contractions.end(); ++iter, ++citer) {
          vector<double>::const_iterator eiter = exponents.begin();
          double denom = 1.0;
          for (int ii = 2; ii <= i; ++ii) denom *= 2 * ii - 1;
          for (vector<double>::iterator diter = iter->begin(); diter != iter->end(); ++diter, ++eiter) {
            *diter *= ::pow(2.0 * *eiter / PI, 0.75) * ::pow(::sqrt(4.0 * *eiter), static_cast<double>(i)) / ::sqrt(denom);
          }
  
          vector<vector<double> > cont(1, *iter);
          vector<pair<int, int> > cran(1, *citer);
          RefShell current(new Shell(spherical_, position_, i, exponents, cont, cran));
          vector<RefShell> cinp(2, current); 
          OverlapBatch coverlap(cinp);
          coverlap.compute();
          const double scal = 1.0 / ::sqrt((coverlap.data())[0]);
          for (vector<double>::iterator diter = iter->begin(); diter != iter->end(); ++diter) *diter *= scal;
        } 

        RefShell currentbatch(new Shell(spherical_, position_, i, exponents, contractions, contranges));
        shells_.push_back(currentbatch);
        lmax_ = i;
      }

    } // end of batch loop 
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


void Atom::print_basis() const {
  for(vector<RefShell>::const_iterator iter = shells_.begin(); iter != shells_.end(); ++iter)
     cout << (*iter)->show() << endl; 
}


void Atom::print() const {
  cout << "  " + name_ << fixed << setprecision(10) <<
      setw(15) << position_[0] <<
      setw(15) << position_[1] <<
      setw(15) << position_[2] <<  endl;
}

