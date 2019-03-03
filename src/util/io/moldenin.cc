//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenin.cc
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

#include <streambuf>
#include <src/util/io/moldenin.h>
#include <src/util/atommap.h>
#include <src/util/constants.h>

#if defined(__GNUC__) && ((__GNUC__ == 4 && __GNUC_MINOR__ >= 9) || __GNUC__ > 4)
#include <regex>
#else
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_search;
using boost::regex_replace;
using boost::cmatch;
#endif

using namespace bagel;
using namespace std;

struct istreambuf : public basic_streambuf<char, std::char_traits<char>> {
    istreambuf(char* buffer, streamsize length) {
        setg(buffer, buffer, buffer + length);
    }
    pos_type seekpos(pos_type sp, ios_base::openmode which) override {
      return seekoff(sp - pos_type(off_type(0)), std::ios_base::beg, which);
    }
    pos_type seekoff(off_type off, ios_base::seekdir dir, ios_base::openmode which = ios_base::in) override {
      if (dir == ios_base::cur)
        gbump(off);
      else if (dir == ios_base::end)
        setg(eback(), egptr() + off, egptr());
      else if (dir == std::ios_base::beg)
        setg(eback(), eback() + off, egptr());
      return gptr() - eback();
    }
};

/************************************************************************************
************************************************************************************/

MoldenIn::MoldenIn(const string filename, const bool is_spherical) : MoldenIO(filename), is_spherical_(is_spherical) {
  compute_transforms();

}

MoldenIn& MoldenIn::operator>> (vector<shared_ptr<const Atom>>& atoms) {
  atoms = atoms_;

  return *this;
}


void MoldenIn::read_mos(MOInfo& moinfo) {
  shared_ptr<const Geometry> geom = moinfo.geom;
  auto coeff = make_shared<Coeff>(geom);
  shared_ptr<Coeff> coeffB;
  if (has_beta()) {
    coeffB = make_shared<Coeff>(geom);
  }

  vector<int> atom_offsets;
  for (auto& ioff : geom->offsets())
    atom_offsets.push_back(!ioff.empty() ? ioff.front(): 0);

  double* idataA = coeff->data();
  double* idataB = coeffB ? coeffB->data() : nullptr;
  size_t offsetA = 0;
  size_t offsetB = 0;
  int n = 0;
  int nA = 0;
  int nB = 0;
  for (auto& imo : mo_coefficients_) {
    double* idata;
    if (mo_spin_[n++] == 0) {
      idata = idataA+offsetA;
      offsetA += geom->nbasis();
      nA++;
    } else {
      idata = idataB+offsetB;
      offsetB += geom->nbasis();
      nB++;
    }
    moread_util(imo.data(), idata, atom_offsets);
  }

  assert(!coeffB || coeff->mdim() == coeffB->mdim());
  moinfo.coeff = make_shared<Coeff>(coeff->slice(0, nA));
  if (coeffB)
    moinfo.coeffB = make_shared<Coeff>(coeffB->slice(0, nB));
}


void MoldenIn::read_mos_complex(MOInfo& moinfo) {
  shared_ptr<const Geometry> geom = moinfo.geom;
  auto coeff = make_shared<ZCoeff>(geom);
  if (has_beta())
    throw runtime_error("Complex Molden interface hasn't been implemented for open-shell cases");

  vector<int> atom_offsets;
  for (auto& ioff : geom->offsets())
    atom_offsets.push_back(!ioff.empty() ? ioff.front(): 0);

  complex<double>* idataA = coeff->data();
  size_t offsetA = 0;
  for (auto& imo : mo_coefficients_complex_) {
    complex<double>* idata = idataA+offsetA;
    offsetA += geom->nbasis();
    moread_util(imo.data(), idata, atom_offsets);
  }

  moinfo.zcoeff = coeff;
}


void MoldenIn::read_mos_relativistic(MOInfo& moinfo) {
  shared_ptr<const Geometry> geom = moinfo.geom;
  const int nbasis = geom->nbasis();
  auto coeff = make_shared<ZCoeff_Striped>(nbasis*4, false, nbasis*2, 0, 0, 0);
  assert(has_beta());

  vector<int> atom_offsets;
  for (auto& ioff : geom->offsets())
    atom_offsets.push_back(!ioff.empty() ? ioff.front(): 0);

  complex<double>* idataA = coeff->data();
  size_t offsetA = 0;
  for (auto& imo : mo_coefficients_relativistic_) {
    vector<molden_impl::complex4> tmp(nbasis);
    moread_util(imo.data(), tmp.data(), atom_offsets);

    complex<double>* idata = idataA+offsetA;
    for (size_t i = 0; i != nbasis; ++i) {
      idata[i]          = tmp[i].data[0];
      idata[i+nbasis]   = tmp[i].data[1];
      idata[i+nbasis*2] = tmp[i].data[2];
      idata[i+nbasis*3] = tmp[i].data[3];
    }
    offsetA += nbasis*4;
  }

  moinfo.relcoeff = coeff;
}


MoldenIn& MoldenIn::operator>> (MOInfo& moinfo) {
  assert(has_mo());
  // read MOs here
  if (!mo_coefficients_.empty())
    read_mos(moinfo);
  else if (!mo_coefficients_complex_.empty())
    read_mos_complex(moinfo);
  else if (!mo_coefficients_relativistic_.empty())
    read_mos_relativistic(moinfo);
  else
    throw logic_error("No MO found; this should not happen");

  shared_ptr<const Geometry> geom = moinfo.geom;

  // output area (to be assigned to the membres of moinfo)
  auto eig   = make_shared<VectorB>(moinfo.norb());
  auto occup = make_shared<VectorB>(moinfo.norb());
  shared_ptr<VectorB> eigB, occupB;
  if (has_beta()) {
    eigB   = eig->clone();
    occupB = occup->clone();
  }

  if (mo_eig_.size() != mo_occup_.size() || mo_eig_.size() != mo_spin_.size())
    throw runtime_error("MoldenIn - inconsistency");

  for (int i = 0, ia = 0; i != mo_spin_.size(); ++i) { 
    if (mo_spin_[i] == 0) {
      (*eig)[ia] = mo_eig_[i];
      (*occup)[ia] = mo_occup_[i];
      ++ia;
    } else {
      assert(!!eigB && !!occupB); 
      (*eigB)[i-ia] = mo_eig_[i];
      (*occupB)[i-ia] = mo_occup_[i];
    }
  }
  assert(!eigB || eig->size() == eigB->size()); 


  moinfo.eig = move(*eig);
  moinfo.occup = move(*occup);
  if (eigB) {
    moinfo.eigB = move(*eigB);
    moinfo.occupB = move(*occupB);
  }
  return *this;
}


void MoldenIn::read() {
  /************************************************************
  *  Set up variables that will contain the organized info    *
  ************************************************************/
  int num_atoms = 0;

  /* Atom positions */
  vector<array<double,3>> positions;
  /* Atom names */
  vector<string> names;
  /* Atom charges */
  vector<double> charges;
  /* Map atom number to basis info */
  map<int, vector<tuple<string, vector<double>, vector<double>>>> basis_info;

  /************************************************************
  *  Set up "global" regular expressions                      *
  ************************************************************/
  regex gto_re("\\[GTO\\]");
  regex atoms_re("\\[[Aa][Tt][Oo][Mm][Ss]\\]");
  regex mo_re("\\[MO\\]");
  regex other_re("\\[\\w+\\]");

  cmatch matches;

  /************************************************************
  *  Booleans to check and make sure each important section   *
  *  was found                                                *
  ************************************************************/
  bool found_atoms = false;
  bool found_gto = false;

  /************************************************************
  * Extra variables                                           *
  ************************************************************/
  double scale;
  string line; // Contains the current line of the file

  /************************************************************
  * Storing the content in a stringstream                     *
  ************************************************************/
  ifstream ifs(filename_);
  if (!ifs.is_open())
    throw runtime_error("Molden input file not found");
  ifs.seekg(0, ios::end);
  const streampos length = ifs.tellg();
  ifs.seekg(0, ios::beg);

  unique_ptr<char[]> buffer(new char[length]);
  ifs.read(buffer.get(), length);
  istreambuf isbuf(buffer.get(), length);
  istream is(&isbuf);

  /************************************************************
  *  An inelegant search to check for the 5D keyword          *
  ************************************************************/
  {
    cartesian_ = true;
    regex _5d_re("\\[5[Dd]\\]");
    regex _5d7f_re("\\[5[Dd]7[Ff]\\]");
    regex _5d7f9g_re("\\[5[Dd]7[Ff]9[Gg]\\]");
    while (!is.eof()) {
      getline(is, line);
      if (regex_search(line,_5d_re) || regex_search(line,_5d7f_re) || regex_search(line,_5d7f9g_re)) {
        cartesian_ = false;
        break;
      }
    }
  }

  /************************************************************
  *  Open input stream                                        *
  ************************************************************/

  is.clear();
  is.seekg(0, is.beg);
  getline(is, line);

  while (!is.eof()) {
    if (regex_search(line,atoms_re)) {
      regex ang_re("Angs");
      regex atoms_line("([a-zA-Z]{1,2})[0-9]*\\s+\\d+\\s+([+-]?[0-9]*[.]?[0-9]+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");

      scale = regex_search(line, ang_re) ? (1.0/au2angstrom__) : 1.0;

      getline(is, line);
      while (!regex_search(line, other_re) && !is.eof()) {
        if (regex_search(line.c_str(), matches, atoms_line)) {
          ++num_atoms;

          const string nm(matches[1].first, matches[1].second);
          names.push_back(nm);

          const string charge_str(matches[2].first, matches[2].second);
          const string x_str(matches[3].first, matches[3].second);
          const string y_str(matches[4].first, matches[4].second);
          const string z_str(matches[5].first, matches[5].second);

          array<double,3> pos;
          pos[0] = lexical_cast<double>(x_str)*scale;
          pos[1] = lexical_cast<double>(y_str)*scale;
          pos[2] = lexical_cast<double>(z_str)*scale;
          positions.push_back(pos);

          charges.push_back(lexical_cast<double>(charge_str));

          getline(is, line);
        }
        else { getline(is, line); }
      }
      found_atoms = true;

    } else if (regex_search(line,gto_re)) {
      getline(is, line);

      regex atom_line("(\\d+)\\s*\\S*");
      regex shell_line("([spdfghij])\\s+(\\d+)\\s*\\S*");
      regex exp_line("(\\S+)\\s+(\\S+)");
      regex Dd("[Dd]");

      while (!regex_search(line, other_re) && !is.eof()) {
        /* This line should be a new atom */
        if (!regex_search(line.c_str(), matches, atom_line)) {
           getline(is, line); continue;
        }
        vector<tuple<string,vector<double>,vector<double>>> atom_basis_info;

        const string atom_no_str(matches[1].first, matches[1].second);
        const int atom_no = lexical_cast<int>(atom_no_str);

        gto_order_.push_back(atom_no);

        vector<int> atomic_shell_order;
        AtomMap atommap;

        getline(is, line);

        while (regex_search(line.c_str(), matches, shell_line)) {
          /* Now it should be a new angular shell */
          const string ang_l(matches[1].first, matches[1].second);
          atomic_shell_order.push_back(atommap.angular_number(ang_l));

          const string num_exp_string(matches[2].first, matches[2].second);
          const int num_exponents = lexical_cast<int>(num_exp_string);

          vector<double> exponents;
          vector<double> coefficients;

          for (int i = 0; i < num_exponents; ++i) {
            getline(is, line);

            regex_search(line.c_str(), matches, exp_line);

            string exp_string(matches[1].first, matches[1].second);
            string coeff_string(matches[2].first, matches[2].second);

            exp_string = regex_replace(exp_string, Dd, "E");
            coeff_string = regex_replace(coeff_string, Dd, "E");

            const double exponent = lexical_cast<double>(exp_string);
            const double coeff = lexical_cast<double>(coeff_string);

            exponents.push_back(exponent);
            coefficients.push_back(coeff);
          }

          /************************************************************
          *  Right now I have evec and cvec which have the newly      *
          *  exponents and coefficients                               *
          ************************************************************/
          atom_basis_info.push_back(make_tuple(ang_l, exponents, coefficients));

          getline(is, line);
        }
        //atom_basis_info to basis_info
        shell_orders_.push_back(atomic_shell_order);
        basis_info.insert(pair<int,vector<tuple<string,vector<double>,vector<double>>>>(atom_no, atom_basis_info));
      }
      found_gto = true;

    } else if (regex_search(line, mo_re)) {
      if (!found_gto) {
        throw runtime_error("MO section found before GTO section. Check Molden file.");
      }
      /* Not used at the moment. Maybe later.
      regex sym_re("Sym=\\s+(\\S+)"); */
      regex ene_re("Ene=\\s*(\\S+)");
      regex spin_re("Spin=\\s*(\\w+)");
      regex occup_re("Occup=\\s*(\\S+)");
      regex coeff_re("\\d+\\s+(\\S+)\\s*$");
      regex coeff_re_comp("\\d+\\s+(\\(\\s*\\S+\\s*,\\s*\\S+\\s*\\))\\s*$");
      regex coeff_re_rel("\\d+\\s+(\\(\\s*\\S+\\s*,\\s*\\S+\\s*\\))\\s*(\\(\\s*\\S+\\s*,\\s*\\S+\\s*\\))\\s*(\\(\\s*\\S+\\s*,\\s*\\S+\\s*\\))\\s*(\\(\\s*\\S+\\s*,\\s*\\S+\\s*\\))\\s*$");

      getline(is, line);
      while (!regex_search(line, other_re) && !is.eof()) {
        if (regex_search(line.c_str(), matches, ene_re)) {
          string estring(matches[1].first, matches[1].second);
          mo_eig_.push_back(lexical_cast<double>(estring));
        }
        if (regex_search(line.c_str(), matches, occup_re)) {
          string estring(matches[1].first, matches[1].second);
          mo_occup_.push_back(lexical_cast<double>(estring));
        }
        if (regex_search(line.c_str(), matches, spin_re)) {
          string estring(matches[1].first, matches[1].second);
          if (estring == "Alpha")
            mo_spin_.push_back(0);
          else if (estring == "Beta")
            mo_spin_.push_back(1);
          else
            throw runtime_error("MoldenIn - spin has to be either Alpha or Beta");
        }
        getline(is, line);
        if (is.eof()) break;

        while (!(regex_search(line.c_str(), coeff_re) || regex_search(line.c_str(), coeff_re_comp) || regex_search(line.c_str(), coeff_re_rel))) {
          if (regex_search(line.c_str(), matches, ene_re)) {
            string estring(matches[1].first, matches[1].second);
            mo_eig_.push_back(lexical_cast<double>(estring));
          }
          if (regex_search(line.c_str(), matches, occup_re)) {
            string estring(matches[1].first, matches[1].second);
            mo_occup_.push_back(lexical_cast<double>(estring));
          }
          if (regex_search(line.c_str(), matches, spin_re)) {
            string estring(matches[1].first, matches[1].second);
            if (estring == "Alpha")
              mo_spin_.push_back(0);
            else if (estring == "Beta")
              mo_spin_.push_back(1);
            else
              throw runtime_error("MoldenIn - spin has to be either Alpha or Beta");
          }
          getline(is, line);
        }

        const bool real = regex_search(line.c_str(), coeff_re) && !regex_search(line.c_str(), coeff_re_comp);
        const bool comp = regex_search(line.c_str(), coeff_re_comp); 
        const bool rel = regex_search(line.c_str(), coeff_re_rel); 
        assert(static_cast<int>(real) + static_cast<int>(comp) + static_cast<int>(rel) == 1);

        if (real) {
          vector<double> movec;
          while (regex_search(line.c_str(), matches, coeff_re)) {
            string mo_string(matches[1].first, matches[1].second);
            movec.push_back(lexical_cast<double>(mo_string));
            getline(is, line);
          }
          mo_coefficients_.push_back(movec);
        } else if (comp) {
          vector<complex<double>> movec;
          while (regex_search(line.c_str(), matches, coeff_re_comp)) {
            string mo_string(matches[1].first, matches[1].second);
            movec.push_back(lexical_cast<complex<double>>(mo_string));
            getline(is, line);
          }
          mo_coefficients_complex_.push_back(movec);
        } else if (rel) {
          vector<molden_impl::complex4> movec;
          while (regex_search(line.c_str(), matches, coeff_re_rel)) {
            molden_impl::complex4 tmp;
            for (int i = 1; i <= 4; ++i) {
              string mo_string(matches[i].first, matches[i].second);
              tmp.data[i-1] = lexical_cast<complex<double>>(mo_string);
            }
            movec.push_back(tmp);
            getline(is, line);
          }
          mo_coefficients_relativistic_.push_back(movec);

        } else {
          throw logic_error("should not happen");
        }
      }
      if (mo_eig_.size() != mo_coefficients_.size() && mo_eig_.size() != mo_coefficients_complex_.size()
       && mo_eig_.size() != mo_coefficients_relativistic_.size())
        throw logic_error("Eigenvalues and eigenvectors are not consistent");
    } else {
      getline(is, line);
    }
  }

  /************************************************************
  *  Check to make sure all the necessary information was     *
  *  found.                                                   *
  ************************************************************/
  if (!(found_atoms && found_gto)) {
     string message("Section not found in Molden file: ");
     if (!found_atoms) { message += "atoms "; }
     if (!found_gto)   { message += "GTO"; }
     throw runtime_error(message);
  }

  vector<shared_ptr<const Atom>> all_atoms;

  /* Assuming the names and positions vectors are in the right order */
  auto iname = names.begin();
  auto piter = positions.begin();
  auto citer = charges.begin();
  for (int i = 0; i < num_atoms; ++i, ++iname, ++piter, ++citer) {
    vector<tuple<string, vector<double>, vector<double>>> binfo = basis_info.find(i+1)->second;
    if (i == num_atoms)
      throw runtime_error("It appears an atom was missing in the GTO section. Check your file");

    /* For each atom, I need to make an atom object and stick it into a vector */
    const string lname = to_lower(*iname);
    shared_ptr<const Atom> at;
    if (lname != "q") {
      at = make_shared<Atom>(is_spherical_, lname, *piter, binfo, "molden");
      if (fabs(at->atom_charge() - *citer) > 1.0e-16)
        throw runtime_error("MoldenIn failed - inconsistent atom charge");
    } else {
      at = make_shared<Atom>(is_spherical_, lname, *piter, *citer);
    }
    all_atoms.push_back(at);
  }

  atoms_ = all_atoms;

}
