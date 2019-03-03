//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: atom.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#include <src/molecule/atom.h>
#include <src/util/math/quatern.h>
#include <src/integral/os/overlapbatch.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

static const AtomMap atommap;

Atom::Atom(shared_ptr<const PTree> inp, const bool spherical, const bool angstrom, const pair<string, shared_ptr<const PTree>> defbas,
           shared_ptr<const PTree> elem, const bool aux, const bool ecp, const bool default_finite)
: spherical_(spherical), use_ecp_basis_(false), basis_(inp->get<string>(!aux ? "basis" : "df_basis", defbas.first)) {
  name_ = to_lower(inp->get<string>("atom"));
  if (basis_.find("ecp") != string::npos) use_ecp_basis_ = true;

  if (elem)
    for (auto& i : *elem) {
      const string key = to_lower(i->key());
      if (name_ == key) basis_ = i->data();
    }
  atom_number_ = atommap.atom_number(name_);

  position_ = inp->get_array<double,3>("xyz");

  for (auto& i : position_) i /= angstrom ? au2angstrom__ : 1.0;

  if (name_ == "q") {
    atom_charge_ = inp->get<double>("charge");
    nbasis_ = 0;
    lmax_ = 0;
  } else {
    shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
    string na = name_;
    na[0] = toupper(na[0]);
    (!use_ecp_basis_) ? basis_init(basisset->get_child(na)) : basis_init_ECP(basisset->get_child(na));
    if (!use_ecp_basis_ && ecp) {
      ecp_parameters_ = make_shared<const ECP>();
      so_parameters_ = make_shared<const SOECP>();
      use_ecp_basis_ = true;
    }
  }
  atom_exponent_ = inp->get<double>("exponent", default_finite ? atommap.nuclear_exponent(name_) : 0.0);

  mass_ = inp->get<double>("mass", atommap.averaged_mass(name_));
}


// constructor that uses the old atom and basis
Atom::Atom(const Atom& old, const bool spherical, const string bas, const pair<string, shared_ptr<const PTree>> defbas, shared_ptr<const PTree> elem)
 : spherical_(spherical), name_(old.name_), position_(old.position_), use_ecp_basis_(old.use_ecp_basis_), ecp_parameters_(old.ecp_parameters_),
   so_parameters_(old.so_parameters_), atom_number_(old.atom_number_), atom_charge_(old.atom_charge_),
   atom_exponent_(old.atom_exponent_), mass_(old.mass_), basis_(bas) {

  if (basis_.find("ecp") != string::npos) use_ecp_basis_ = true;
  if (name_ == "q") {
    nbasis_ = 0;
    lmax_ = 0;
  } else {
    if (elem)
      for (auto& i : *elem) {
        const string key = to_lower(i->key());
        if (name_ == key) basis_ = i->data();
      }
    string na = name_;
    shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
    na[0] = toupper(na[0]);
    (!use_ecp_basis_) ? basis_init(basisset->get_child(na)) : basis_init_ECP(basisset->get_child(na));
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

void Atom::basis_init_ECP(shared_ptr<const PTree> basis) {

  for (auto& ibas : *basis) {
    try
    {
      basis_init(ibas->get_child("valence"));
    }
    catch (const exception &err)
    {
      cout << err.what() << endl;
      throw runtime_error("ECP basis set file has the wrong format!");
    }
    const int ncore = ibas->get<int>("ncore");
    const shared_ptr<const PTree> core = ibas->get_child("core");
    vector<tuple<string, vector<double>, vector<double>, vector<int>>> basis_info;

    for (auto& ibcore : *core) {

      const string ang = ibcore->get<string>("ecp_angular");
      const shared_ptr<const PTree> exp = ibcore->get_child("ecp_exp");
      vector<double> exponents;

      for (auto& p : *exp)
        exponents.push_back(lexical_cast<double>(p->data()));

      const shared_ptr<const PTree> coef = ibcore->get_child("ecp_coef");
      vector<double> coefficients;

      for (auto& c : *coef)
          coefficients.push_back(lexical_cast<double>(c->data()));

      const shared_ptr<const PTree> r_p = ibcore->get_child("ecp_r");
      vector<int> r_power;

      for (auto& r : *r_p)
          r_power.push_back(lexical_cast<int>(r->data()));

      basis_info.push_back(make_tuple(ang, exponents, coefficients, r_power));
    }

    construct_shells_ECP(ncore, basis_info);


    const shared_ptr<const PTree> soecp = ibas->get_child_optional("so");
    if (soecp) {
      const shared_ptr<const PTree> soecp = ibas->get_child("so");
      vector<tuple<string, vector<double>, vector<double>, vector<int>>> sobasis_info;
      for (auto& ibcore : *soecp) {

        const string ang = ibcore->get<string>("ecp_angular");
        const shared_ptr<const PTree> exp = ibcore->get_child("ecp_exp");
        vector<double> exponents;

        for (auto& p : *exp)
          exponents.push_back(lexical_cast<double>(p->data()));

        const shared_ptr<const PTree> coef = ibcore->get_child("ecp_coef");
        vector<double> coefficients;

        for (auto& c : *coef)
            coefficients.push_back(lexical_cast<double>(c->data()));

        const shared_ptr<const PTree> r_p = ibcore->get_child("ecp_r");
        vector<int> r_power;

        for (auto& r : *r_p)
            r_power.push_back(lexical_cast<int>(r->data()));

        sobasis_info.push_back(make_tuple(ang, exponents, coefficients, r_power));
      }

      construct_shells_SOECP(sobasis_info);
    } else {
      so_parameters_ = make_shared<const SOECP>();
    }
  }

}

Atom::Atom(const Atom& old, const array<double, 3>& displacement)
: spherical_(old.spherical_), name_(old.name()), use_ecp_basis_(old.use_ecp_basis_), ecp_parameters_(old.ecp_parameters_),
  so_parameters_(old.so_parameters_),  atom_number_(old.atom_number()), atom_charge_(old.atom_charge()),
  atom_exponent_(old.atom_exponent()), mass_(old.mass_), nbasis_(old.nbasis()), lmax_(old.lmax()), basis_(old.basis_) {

  assert(displacement.size() == 3 && old.position().size() == 3);
  const array<double,3> opos = old.position();
  position_ = array<double,3>{{displacement[0]+opos[0], displacement[1]+opos[1], displacement[2]+opos[2]}};

  const vector<shared_ptr<const Shell>> old_shells = old.shells();
  for(auto& s : old_shells)
    shells_.push_back(s->move_atom(displacement));
}

Atom::Atom(const string nm, const string bas, const vector<shared_ptr<const Shell>> shell,
                                              const vector<shared_ptr<const Shell_ECP>> shell_ECP, const int ncore, const int maxl)
: name_(nm), shells_(shell), use_ecp_basis_(true), ecp_parameters_(make_shared<const ECP>(ncore, maxl, shell_ECP)), atom_number_(atommap.atom_number(nm)), basis_(bas) {
  spherical_ = shells_.front()->spherical();
  position_ = shells_.front()->position();

  common_init();
  atom_exponent_ = 0.0;
  mass_ = atommap.averaged_mass(name_);
}

Atom::Atom(const string nm, const string bas, const vector<shared_ptr<const Shell>> shell, const shared_ptr<const ECP> ecp_param)
: name_(nm), shells_(shell), use_ecp_basis_(true), ecp_parameters_(ecp_param), atom_number_(atommap.atom_number(nm)), basis_(bas) {
  spherical_ = shells_.front()->spherical();
  position_ = shells_.front()->position();

  common_init();
  atom_exponent_ = 0.0;
  mass_ = atommap.averaged_mass(name_);
}


Atom::Atom(const string nm, const string bas, vector<shared_ptr<const Shell>> shell)
: name_(nm), shells_(shell), use_ecp_basis_(false), atom_number_(atommap.atom_number(nm)), basis_(bas) {
  spherical_ = shells_.front()->spherical();
  position_ = shells_.front()->position();

  common_init();
  atom_exponent_ = 0.0;
  mass_ = atommap.averaged_mass(name_);
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, const string bas, const pair<string, shared_ptr<const PTree>> defbas, shared_ptr<const PTree> elem)
 : spherical_(sph), name_(nm), position_(p), use_ecp_basis_(false), atom_number_(atommap.atom_number(nm)), basis_(bas) {

  if (elem)
    for (auto& i : *elem) {
      const string key = to_lower(i->key());
      if (name_ == key) basis_ = i->data();
    }

  string na = name_;
  na[0] = toupper(na[0]);
  shared_ptr<const PTree> basisset = (basis_ == defbas.first) ? defbas.second : PTree::read_basis(basis_);
  if (basis_.find("ecp") != string::npos) use_ecp_basis_ = true;
  (!use_ecp_basis_) ? basis_init(basisset->get_child(na)) : basis_init_ECP(basisset->get_child(na));

  atom_exponent_ = 0.0;
  mass_ = atommap.averaged_mass(name_);
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, vector<tuple<string, vector<double>, vector<double>>> in, const string bas)
 : spherical_(sph), name_(nm), position_(p), use_ecp_basis_(false), atom_number_(atommap.atom_number(nm)), basis_(bas) {

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
  mass_ = atommap.averaged_mass(name_);
}


Atom::Atom(const bool sph, const string nm, const array<double,3>& p, const double charge)
: spherical_(sph), name_(nm), position_(p), use_ecp_basis_(false), atom_number_(atommap.atom_number(nm)), atom_charge_(charge), nbasis_(0), lmax_(0), basis_("") {
  atom_exponent_ = 0.0;
  mass_ = atommap.averaged_mass(name_);
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

  for (int i = 0; i <= atommap.max_angular_number(); ++i) {
    vector<vector<double>> contractions;
    vector<pair<int, int>> contranges;
    vector<double> exponents;

    int offset = 0;
    for (auto biter = in.begin(); biter != in.end(); ++biter) {

      // check the angular number
      if (atommap.angular_number(get<0>(*biter)) != i) continue;

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
        contranges.push_back({offset + zerostart, offset + current.size() - zeroend});
        assert(offset + zerostart <= offset + current.size() - zeroend);
      }
      const vector<double> exp = get<1>(*biter);
      if (biter+1 == in.end() || exp != get<1>(*(biter+1)) || atommap.angular_number(get<0>(*(biter+1))) != i) {
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
          *diter *= pow(2.0 * *eiter / pi__, 0.75) * pow(sqrt(4.0 * *eiter), static_cast<double>(i)) / sqrt(denom);

        if (basis_ != "molden") {
          vector<vector<double>> cont {*iter};
          vector<pair<int, int>> cran {*citer};
          auto current = make_shared<const Shell>(spherical_, position_, i, exponents, cont, cran);
          array<shared_ptr<const Shell>,2> cinp {{ current, current }};
          OverlapBatch coverlap(cinp);
          coverlap.compute();
          const double scal = 1.0 / sqrt((coverlap.data())[0]);
          for (auto& d : *iter) d *= scal;
        }
      }

      shells_.push_back(make_shared<Shell>(spherical_, position_, i, exponents, contractions, contranges));
      lmax_ = i;
    }

  } // end of batch loop

  // TODO size is not optimized!
  split_shells(40);

}

void Atom::construct_shells_ECP(const int ncore, vector<tuple<string, vector<double>, vector<double>, vector<int>>> in) {

  vector<shared_ptr<const Shell_ECP>> shells_ECP;

  int maxl = 0;
  for (auto& biter : in) {
    const int l = atommap.angular_number(get<0>(biter));
    if (l > maxl) maxl = l;
    vector<double> exponents = get<1>(biter);
    const vector<double> coefficients = get<2>(biter);
    vector<int> r_power = get<3>(biter);

    vector<double> new_coefficients(coefficients);
    int nzero = 0;
    for (auto ic = coefficients.begin(); ic != coefficients.end(); ++ic)
       if (*ic == 0.0) {
         auto pos = std::distance(coefficients.begin(), ic) - nzero;
         exponents.erase(exponents.begin()+pos);
         r_power.erase(r_power.begin()+pos);
         new_coefficients.erase(new_coefficients.begin()+pos);
         ++nzero;
       }
    assert(exponents.size() == r_power.size() && r_power.size() == new_coefficients.size());

    if (!exponents.empty()) shells_ECP.push_back(make_shared<const Shell_ECP>(position_, l , exponents, new_coefficients, r_power));

  }

  ecp_parameters_ = (shells_ECP.empty()) ? make_shared<const ECP>() : make_shared<const ECP>(ncore, maxl, shells_ECP);

}

void Atom::construct_shells_SOECP(vector<tuple<string, vector<double>, vector<double>, vector<int>>> in) {

  vector<shared_ptr<const Shell_ECP>> shells_SO;

  for (auto& biter : in) {
    const int l = atommap.angular_number(get<0>(biter));
    vector<double> exponents = get<1>(biter);
    const vector<double> coefficients = get<2>(biter);
    vector<int> r_power = get<3>(biter);

    vector<double> new_coefficients(coefficients);
    int nzero = 0;
    for (auto ic = coefficients.begin(); ic != coefficients.end(); ++ic)
       if (*ic == 0.0) {
         auto pos = std::distance(coefficients.begin(), ic) - nzero;
         exponents.erase(exponents.begin()+pos);
         r_power.erase(r_power.begin()+pos);
         new_coefficients.erase(new_coefficients.begin()+pos);
         ++nzero;
       }
    assert(exponents.size() == r_power.size() && r_power.size() == new_coefficients.size());

    if (!exponents.empty()) shells_SO.push_back(make_shared<const Shell_ECP>(position_, l , exponents, new_coefficients, r_power));
  }

  so_parameters_ = (shells_SO.empty()) ? make_shared<const SOECP>() : make_shared<const SOECP>(shells_SO);

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
  for (auto& i : shells_)
    cout << i->show() << endl;

  if (ecp_parameters_)
    ecp_parameters_->print();
}


void Atom::print() const {
  string tmp = name_;
  tmp[0] = toupper(tmp[0]);
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
  return atan2(rot.norm(), ap.dot_product(bp)) * rad2deg__;
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
  return atan2(b2.norm()*b1.dot_product(b23), b12.dot_product(b23)) * rad2deg__;
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
  atom->reset_shells(rshells);
  return atom;
}


shared_ptr<const Atom> Atom::relativistic(const array<double,3>& magnetic_field, bool london) const {
  // basically the same
  // except for shells_
  vector<shared_ptr<const Shell>> rshells;
  for (auto& i : shells_) {
    auto tmp = make_shared<Shell>(*i);
    tmp->init_relativistic(magnetic_field, london);
    rshells.push_back(tmp);
  }
  auto atom = make_shared<Atom>(*this);
  atom->shells_ = rshells;
  return atom;
}


shared_ptr<const Atom> Atom::apply_magnetic_field(const array<double,3>& magnetic_field, const bool london) const {

  auto atom = make_shared<Atom>(*this);
  if (london) {
    atom->vector_potential_[0] = 0.5*(magnetic_field[1]*position_[2] - magnetic_field[2]*position_[1]);
    atom->vector_potential_[1] = 0.5*(magnetic_field[2]*position_[0] - magnetic_field[0]*position_[2]);
    atom->vector_potential_[2] = 0.5*(magnetic_field[0]*position_[1] - magnetic_field[1]*position_[0]);
  } else {
    atom->vector_potential_ = array<double,3>{{0.0, 0.0, 0.0}};
  }

  // basically the same
  // except for shells_
  vector<shared_ptr<const Shell>> mshells;
  for (auto& i : shells_) {
    auto tmp = make_shared<Shell>(*i);
    tmp->add_phase(atom->vector_potential_, magnetic_field, london);
    mshells.push_back(tmp);
  }
  atom->shells_ = mshells;

  return atom;
}


double Atom::radius() const { return atommap.radius(name_); }
double Atom::cov_radius() const { return atommap.cov_radius(name_); }


shared_ptr<const Atom> Atom::uncontract() const {
  auto atom = make_shared<Atom>(*this);

  vector<shared_ptr<const Shell>> uncshells;
  for (auto& i : shells_)
    uncshells.push_back(i->uncontract());
  atom->reset_shells(uncshells);

  return atom;
}

void Atom::reset_shells(vector<shared_ptr<const Shell>> rshells) {
  shells_ = rshells;
  common_init();
}
