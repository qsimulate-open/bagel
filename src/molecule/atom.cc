//
// BAGEL - Parallel electron correlation program.
// Filename: atom.cc
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


#include <src/molecule/atom.h>
#include <src/math/quatern.h>
#include <src/integral/os/overlapbatch.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

static const AtomMap atommap_;

Atom::Atom(shared_ptr<const PTree> inp, const bool spherical, const bool angstrom, const pair<string, shared_ptr<const PTree>> defbas, std::shared_ptr<const PTree> elem, const bool aux)
: spherical_(spherical), basis_(inp->get<string>(!aux ? "basis" : "df_basis", defbas.first)) {
  name_ = inp->get<string>("atom");
  transform(name_.begin(), name_.end(), name_.begin(), ::tolower);

  if (elem)
    for (auto& i : *elem) {
      string key = i->key();
      transform(key.begin(), key.end(), key.begin(), ::tolower);
      if (name_ == key) basis_ = i->data();
    }
  atom_number_ = atommap_.atom_number(name_);

  position_ = inp->get_array<double,3>("xyz");
  for (auto& i : position_) i *= angstrom ? ang2bohr__ : 1.0;

  if (name_ == "q") {
    atom_charge_ = inp->get<double>("charge");
    nbasis_ = 0;
    lmax_ = 0;
  } else {
    shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
    string na = name_;
    na[0] = toupper(na[0]);
    basis_init(basisset->get_child(na));
  }
  atom_exponent_ = inp->get<double>("exponent", 0.0);
}


// constructor that uses the old atom and basis
Atom::Atom(const Atom& old, const bool spherical, const string bas, const pair<string, shared_ptr<const PTree>> defbas, std::shared_ptr<const PTree> elem)
 : spherical_(spherical), name_(old.name_), position_(old.position_), atom_number_(old.atom_number_), atom_charge_(old.atom_charge_), atom_exponent_(old.atom_exponent_), basis_(bas) {
  if (name_ == "q") {
    nbasis_ = 0;
    lmax_ = 0;
  } else {
    if (elem)
      for (auto& i : *elem) {
        string key = i->key();
        transform(key.begin(), key.end(), key.begin(), ::tolower);
        if (name_ == key) basis_ = i->data();
      }
    string na = name_;
    shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
    na[0] = toupper(na[0]);
    basis_init(basisset->get_child(na));
  }
}


void Atom::basis_init(shared_ptr<const PTree> basis) {
  // basis_info will be used in the construction of Basis_batch
  vector<tuple<string, vector<double>, vector<vector<double>>>> basis_info;

  for (auto& ibas : *basis) {

    const string ang = ibas->get<string>("angular");
    const shared_ptr<const PTree> prim = ibas->get_child("prim");
    vector<double> exponents;

    for (auto& p : *prim)
      exponents.push_back(lexical_cast<double>(p->data()));

    const shared_ptr<const PTree> cont = ibas->get_child("cont");
    vector<vector<double>> coeff;

    for (auto& c : *cont) {
      vector<double> tmp;
      for (auto& cc : *c)
        tmp.push_back(lexical_cast<double>(cc->data()));
      coeff.push_back(tmp);
    }
    basis_info.push_back(make_tuple(ang, exponents, coeff));
  }
  construct_shells(basis_info);
  common_init();
}


Atom::Atom(const Atom& old, const array<double, 3>& displacement)
: spherical_(old.spherical_), name_(old.name()), atom_number_(old.atom_number()), atom_charge_(old.atom_charge()), atom_exponent_(old.atom_exponent()),
  nbasis_(old.nbasis()), lmax_(old.lmax()), basis_(old.basis_) {

  assert(displacement.size() == 3 && old.position().size() == 3);
  const array<double,3> opos = old.position();
  position_ = array<double,3>{{displacement[0]+opos[0], displacement[1]+opos[1], displacement[2]+opos[2]}};

  const vector<shared_ptr<const Shell>> old_shells = old.shells();
  for(auto& s : old_shells)
    shells_.push_back(s->move_atom(displacement));
}


Atom::Atom(const string nm, const string bas, vector<shared_ptr<const Shell>> shell)
: name_(nm), shells_(shell), atom_number_(atommap_.atom_number(nm)), basis_(bas) {
  spherical_ = shells_.front()->spherical();
  position_ = shells_.front()->position();

  common_init();
  atom_exponent_ = 0.0;
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, const string bas, const std::pair<std::string, std::shared_ptr<const PTree>> defbas, std::shared_ptr<const PTree> elem)
 : spherical_(sph), name_(nm), position_(p), atom_number_(atommap_.atom_number(nm)), basis_(bas) {

  if (elem)
    for (auto& i : *elem) {
      string key = i->key();
      transform(key.begin(), key.end(), key.begin(), ::tolower);
      if (name_ == key) basis_ = i->data();
    }
  string na = name_;
  na[0] = toupper(na[0]);
  shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
  basis_init(basisset->get_child(na));

  atom_exponent_ = 0.0;
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, vector<tuple<string, vector<double>, vector<double>>> in)
 : spherical_(sph), name_(nm), position_(p), atom_number_(atommap_.atom_number(nm)), basis_("custom_basis") {

  // tuple
  vector<tuple<string, vector<double>, vector<vector<double>>>> basis_info;
  for (auto& iele : in) {
    vector<double> tmp = get<2>(iele);
    vector<vector<double>> tmp2 = {tmp};
    basis_info.push_back(make_tuple(get<0>(iele), get<1>(iele), tmp2));
  }

  construct_shells(basis_info);
  common_init();
  atom_exponent_ = 0.0;
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, const double charge)
: spherical_(sph), name_(nm), position_(p), atom_number_(atommap_.atom_number(nm)), atom_charge_(charge), nbasis_(0), lmax_(0), basis_("") {
  atom_exponent_ = 0.0;
}



void Atom::common_init() {
  // counting the number of basis functions belonging to this atom
  nbasis_ = accumulate(shells_.begin(), shells_.end(), 0, [](const int& i, const shared_ptr<const Shell>& j) { return i+j->nbasis(); });
  atom_charge_ = static_cast<double>(atom_number_);
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
void Atom::construct_shells(vector<tuple<string, vector<double>, vector<vector<double>>>> in) {

  for (int i = 0; i <= atommap_.max_angular_number(); ++i) {
    vector<vector<double>> contractions;
    vector<pair<int, int>> contranges;
    vector<double> exponents;

    int offset = 0;
    for (auto biter = in.begin(); biter != in.end(); ++biter) {

      // check the angular number
      if (atommap_.angular_number(get<0>(*biter)) != i) continue;

      // contraction coefficient matrix
      const vector<vector<double>> conts = get<2>(*biter);

      // loop over contraction coefficients
      for (auto& current : conts) {

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
      if (biter+1 == in.end() || exp != get<1>(*(biter+1)) || atommap_.angular_number(get<0>(*(biter+1))) != i) {
        exponents.insert(exponents.end(), exp.begin(), exp.end());
        offset += exp.size();
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
          *diter *= ::pow(2.0 * *eiter / pi__, 0.75) * ::pow(::sqrt(4.0 * *eiter), static_cast<double>(i)) / ::sqrt(denom);

        vector<vector<double>> cont(1, *iter);
        vector<pair<int, int>> cran(1, *citer);
        auto current = make_shared<const Shell>(spherical_, position_, i, exponents, cont, cran);
        array<shared_ptr<const Shell>,2> cinp {{ current, current }};
        OverlapBatch coverlap(cinp);
        coverlap.compute();
        const double scal = 1.0 / ::sqrt((coverlap.data())[0]);
        for (auto& d : *iter) d *= scal;
      }

      shells_.push_back(make_shared<const Shell>(spherical_, position_, i, exponents, contractions, contranges));
      lmax_ = i;
    }

  } // end of batch loop

  // shuffle, but deterministic
  // FIXME this breaks the atomic density guess, since it relies on the shell ordering
#if 0
  srand(0);
  random_shuffle(shells_.begin(), shells_.end(), [](const int& i) { return rand()%i; });
#endif

  // TODO size is not optimized!
  split_shells(40);

}


void Atom::split_shells(const size_t batchsize) {
  vector<shared_ptr<const Shell>> out;
  for (auto& i : shells_) {
    const int nbasis = i->nbasis();
    if (nbasis >= batchsize) {
      vector<shared_ptr<const Shell>> tmp = i->split_if_possible(batchsize);
      out.insert(out.end(), tmp.begin(), tmp.end());
    } else {
      out.push_back(i);
    }
  }
  shells_ = out;
}


void Atom::print_basis() const {
  for (auto& i : shells_) cout << i->show() << endl;
}


void Atom::print() const {
  string tmp = name_;
  tmp[0] = ::toupper(tmp[0]);
  cout << "  { \"atom\" : \"" << tmp << "\", \"xyz\" : [" << fixed << setprecision(6) <<
      setw(14) << position_[0] << "," <<
      setw(14) << position_[1] << "," <<
      setw(14) << position_[2] << " ]";
  if (dummy()) {
    cout << ", \"charge\" : " << setw(14) << atom_charge_;
  }
  cout <<  " }," << endl;
}


bool Atom::operator==(const Atom& o) const {
  bool out = true;
  out &= spherical_ == o.spherical_;
  out &= name_ == o.name_;
  out &= fabs(position_[0]-o.position_[0]) < numerical_zero__;
  out &= fabs(position_[1]-o.position_[1]) < numerical_zero__;
  out &= fabs(position_[2]-o.position_[2]) < numerical_zero__;
  out &= shells_.size() == o.shells_.size();

  for (auto i = shells_.begin(), j = o.shells_.begin(); i != shells_.end(); ++i, ++j) out &= (**i) == (**j);

  out &= atom_number_ == o.atom_number_;
  out &= fabs(atom_charge_ - o.atom_charge_) < numerical_zero__;
  out &= nbasis_ == o.nbasis_;
  out &= lmax_ == o.lmax_;

  return out;
}


double Atom::distance(const array<double,3>& o) const {
  double out = 0.0;
  for (int i = 0; i != 3; ++i)
    out += pow(position_[i] - o[i], 2);
  return sqrt(out);
}


array<double,3> Atom::displ(const shared_ptr<const Atom> o) const {
  return array<double,3>{{ o->position_[0]-position_[0], o->position_[1]-position_[1], o->position_[2]-position_[2] }};
}


double Atom::angle(const shared_ptr<const Atom> a, const shared_ptr<const Atom> b) const {
  Quatern<double> ap = a->position();
  Quatern<double> bp = b->position();
  Quatern<double> op = this->position();
  ap -= op;
  bp -= op;
  Quatern<double> rot = ap * bp;
  rot[0] = 0;
  return ::atan2(rot.norm(), ap.dot_product(bp)) * rad2deg__;
}


// Dihedral angle of A-this-O-B
double Atom::dihedral_angle(const shared_ptr<const Atom> a, const shared_ptr<const Atom> o, const shared_ptr<const Atom> b) const {
  Quatern<double> ap = a->position();
  Quatern<double> tp = this->position();
  Quatern<double> op = o->position();
  Quatern<double> bp = b->position();
  // following Wikipedia..
  Quatern<double> b1 = tp - ap;
  Quatern<double> b2 = op - tp;
  Quatern<double> b3 = bp - op;
  Quatern<double> b12 = b1 * b2; b12[0] = 0.0;
  Quatern<double> b23 = b2 * b3; b23[0] = 0.0;
  return ::atan2(b2.norm()*b1.dot_product(b23), b12.dot_product(b23)) * rad2deg__;
}


shared_ptr<const Atom> Atom::relativistic() const {
  // basically the same
  // except for shells_
  vector<shared_ptr<const Shell>> rshells;
  for (auto& i : shells_) {
    auto tmp = make_shared<Shell>(*i);
    tmp->init_relativistic();
    rshells.push_back(tmp);
  }
  auto atom = make_shared<Atom>(*this);
  atom->shells_ = rshells;
  return atom;
}

double Atom::radius() const { return atommap_.radius(name_); }
double Atom::cov_radius() const { return atommap_.cov_radius(name_); }
