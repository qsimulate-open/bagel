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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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
#include <chrono>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <src/scf/geometry.h>
#include <src/scf/atommap.h>
#include <src/molden/molden.h>
#include <src/util/constants.h>
#include <src/util/quatern.h>
#include <src/rysint/libint.h>

using namespace std;
using namespace bagel;

typedef shared_ptr<const Atom> RefAtom;

const static AtomMap atommap_;

Geometry::Geometry(const multimap<string, string> geominfo)
  : spherical_(true), input_(""), lmax_(0) {

  schwarz_thresh_ = read_input<double>(geominfo, "schwarz_thresh", 1.0e-12);
  overlap_thresh_ = read_input<double>(geominfo, "thresh_overlap", 1.0e-8);

  // symmetry
  symmetry_ = read_input<string>(geominfo, "symmetry", "c1");

  // cartesian or not.
  const bool cart = read_input<bool>(geominfo, "cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }

  /* Set up atoms_ */
  basisfile_ = read_input<string>(geominfo, "basis", "");
  if (basisfile_ == "") throw runtime_error("There is no basis specification");
  else if (basisfile_ == "molden") {
    Molden molden(spherical_);
    string molden_file = read_input<string>(geominfo, "molden_in", "");
    if (molden_file == "") throw runtime_error("No molden_in file provided");
    atoms_ = molden.read_geo(molden_file);
  }
  else {
    const bool angstrom = read_input<bool>(geominfo, "angstrom", false);

    pair<multimap<string,string>::const_iterator, multimap<string,string>::const_iterator> bound = geominfo.equal_range("atom");
    for (auto iter = bound.first; iter != bound.second; ++iter) {
      boost::smatch what;
      const boost::regex atom_reg("\\(\\s*([A-Za-z]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+)\\s*");
      auto start = iter->second.begin();
      auto end = iter->second.end();
      if (regex_search(start, end, what, atom_reg)) {
        const string aname(what[1].first, what[1].second);
        const string x_str(what[2].first, what[2].second);
        const string y_str(what[3].first, what[3].second);
        const string z_str(what[4].first, what[4].second);
        const double prefac = angstrom ? ang2bohr__ : 1.0 ;
        array<double,3> positions
          = {{boost::lexical_cast<double>(x_str)*prefac, boost::lexical_cast<double>(y_str)*prefac, boost::lexical_cast<double>(z_str)*prefac}};
        if (aname != "q") {
          // standard atom construction
          auto next = what[0].second;
          const boost::regex end_reg("\\)");
          if (!regex_search(next, end, what, end_reg))
            throw runtime_error("One of the atom lines is corrupt");
          RefAtom catom(new Atom(spherical_, aname, positions, basisfile_));
          atoms_.push_back(catom);
        } else {
          // dummy atom construction
          if (symmetry_ != "c1")
            throw runtime_error("External point charges are only allowed in C1 calculations so far.");
          auto next = what[0].second;
          const boost::regex charge_reg(",\\s*([0-9\\.-]+)\\s*\\)");
          if (regex_search(next, end, what, charge_reg)) {
            const string charge_str(what[1].first, what[1].second);
            const double charge = boost::lexical_cast<double>(charge_str);
            RefAtom catom(new Atom(spherical_, aname, positions, charge));
            atoms_.push_back(catom);
          } else {
            throw runtime_error("No charge specified for dummy atom(s)");
          }
        }
      } else {
        throw runtime_error("One of the atom lines is corrupt");
      }
    }
  }
  if(atoms_.empty()) throw runtime_error("No atoms specified at all");

  /* Set up aux_atoms_ */
  auxfile_ = read_input<string>(geominfo, "df_basis", "");
  if (!auxfile_.empty()) {
    for(auto iatom = atoms_.begin(); iatom != atoms_.end(); ++iatom) {
      if (!(*iatom)->dummy()) {
        RefAtom catom(new Atom(spherical_, (*iatom)->name(), (*iatom)->position(), auxfile_));
        aux_atoms_.push_back(catom);
      } else {
        // we need a dummy atom here to be consistent in gradient computations
        aux_atoms_.push_back(*iatom);
      }
    }
  }

  // Misc
  aux_merged_ = false;

  common_init1();

  print_atoms();

  common_init2(true, overlap_thresh_);

  // static external field
  external_ = vector<double>(3);
  external_[0] = read_input<double>(geominfo, "ex", 0.0);
  external_[1] = read_input<double>(geominfo, "ey", 0.0);
  external_[2] = read_input<double>(geominfo, "ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
}


void Geometry::common_init1() {
  lmax_ = 0;
  aux_lmax_ = 0;
  nbasis_ = 0;
  naux_ = 0;
  nele_ = 0;
  nfrc_ = 0;

  for (auto catom = atoms_.begin(); catom != atoms_.end(); ++catom) {
    nele_ += atommap_.atom_number((*catom)->name());

    int cc = 0;
    vector<int> coffsets;
    for (auto iter = (*catom)->shells().begin(); iter != (*catom)->shells().end(); ++iter) {
      coffsets.push_back(nbasis_ + cc);
      const int ang = (*iter)->angular_number();
      const int angsize = spherical_ ? (2*ang+1) : (ang+1)*(ang+2)/2;
      cc += angsize * (*iter)->num_contracted();
    }
    lmax_ = max(lmax_, (*catom)->lmax());
    nbasis_ += (*catom)->nbasis();
    offsets_.push_back(coffsets);
  }

  if (!auxfile_.empty()) {
    for (auto catom = aux_atoms_.begin(); catom != aux_atoms_.end(); ++catom) {
      int cc = 0;
      vector<int> coffsets;
      for (auto iter = (*catom)->shells().begin(); iter != (*catom)->shells().end(); ++iter) {
        coffsets.push_back(naux_ + cc);
        const int ang = (*iter)->angular_number();
        const int angsize = spherical_ ? (2*ang+1) : (ang+1)*(ang+2)/2;
        cc += angsize * (*iter)->num_contracted();
      }
      aux_lmax_ = max(lmax_, (*catom)->lmax());
      naux_ += (*catom)->nbasis();
      aux_offsets_.push_back(coffsets);
    }
  }
}


void Geometry::common_init2(const bool print, const double thresh, const bool nodf) {
  // symmetry set-up
  plist_ = std::shared_ptr<Petite>(new Petite(atoms_, symmetry_));
  nirrep_ = plist_->nirrep();

  if (print) {
    cout << endl;
    cout << "  Number of basis functions: " << setw(8) << nbasis() << endl;
    cout << "  Number of electrons      : " << setw(8) << nele() << endl << endl;
  }

  nuclear_repulsion_ = compute_nuclear_repulsion();

  if (!auxfile_.empty() && !nodf) {
    if (print) cout << "  Number of auxiliary basis functions: " << setw(8) << naux() << endl << endl;
    cout << "  Since a DF basis is specified, we compute 2- and 3-index integrals:" << endl;
    cout << "    o Being stored without compression. Storage requirement is "
         << setprecision(3) << static_cast<size_t>(naux_)*nbasis()*nbasis()*8.e-9 << " GB" << endl;
    auto tp1 = chrono::high_resolution_clock::now();
    df_ = form_fit<ERIFit>(thresh, true); // true means we construct J^-1/2

    auto tp2 = chrono::high_resolution_clock::now();
    cout << "        elapsed time:  " << setw(10) << setprecision(2) << chrono::duration_cast<chrono::milliseconds>(tp2-tp1).count()*0.001 <<
            " sec." << endl << endl;
  }

}



// suitable for geometry updates in optimization
Geometry::Geometry(const Geometry& o, const vector<double> displ, const multimap<string, string> geominfo, const bool rotate, const bool nodf)
  : spherical_(o.spherical_), input_(o.input_), aux_merged_(o.aux_merged_), basisfile_(o.basisfile_),
    auxfile_(o.auxfile_), symmetry_(o.symmetry_), schwarz_thresh_(o.schwarz_thresh_), external_(o.external_), gamma_(o.gamma_) {

  // first construct atoms using displacements
  auto disp = displ.begin();
  for (auto i = o.atoms_.begin(), j = o.aux_atoms_.begin(); i != o.atoms_.end(); ++i, ++j, disp += 3) {
    array<double,3> cdispl = {{*disp, *(disp+1), *(disp+2)}};
    atoms_.push_back(shared_ptr<Atom>(new Atom(**i, cdispl)));
    aux_atoms_.push_back(shared_ptr<Atom>(new Atom(**j, cdispl)));
  }

  // second find the unique frame.
  // (i) center of chages
  if (rotate) {
    Quatern<double> oc = o.charge_center();
    Quatern<double> mc = charge_center();
    // (2) direction of the first atom
    Quatern<double> oa = o.atoms().front()->position();
    Quatern<double> ma =   atoms().front()->position();
    Quatern<double> od = oa - oc;
    Quatern<double> md = ma - mc;
    // Quaternion that maps md to od.
    od.normalize();
    md.normalize();
    Quatern<double> op = md * od;
    op[0] = 1.0 - op[0];
    op.normalize();
    Quatern<double> opd = op.dagger();

    // first subtract mc, rotate, and then add oc
    vector<shared_ptr<const Atom> > newatoms;
    vector<shared_ptr<const Atom> > newauxatoms;
    for (auto i = atoms_.begin(), j = aux_atoms_.begin(); i != atoms_.end(); ++i, ++j) {
      assert((*i)->position() == (*j)->position());
      Quatern<double> source = (*i)->position();
      Quatern<double> target = op * (source - mc) * opd + oc;
      array<double,3> cdispl = (target - source).ijk();

      newatoms.push_back(shared_ptr<Atom>(new Atom(**i, cdispl)));
      newauxatoms.push_back(shared_ptr<Atom>(new Atom(**j, cdispl)));
    }
    atoms_ = newatoms;
    aux_atoms_ = newauxatoms;
  }

  common_init1();
  overlap_thresh_ = read_input<double>(geominfo, "thresh_overlap", 1.0e-8);
  common_init2(false, overlap_thresh_, nodf);
}


Geometry::Geometry(const Geometry& o, const array<double,3> displ)
  : spherical_(o.spherical_), input_(o.input_), aux_merged_(o.aux_merged_), basisfile_(o.basisfile_),
    auxfile_(o.auxfile_), symmetry_(o.symmetry_), schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_),
    external_(o.external_), gamma_(o.gamma_) {

  // first construct atoms using displacements

  for (auto i = o.atoms_.begin(); i != o.atoms_.end(); ++i) {
    atoms_.push_back(shared_ptr<Atom>(new Atom(**i, displ)));
  }

  for (auto j = o.aux_atoms_.begin(); j != o.aux_atoms_.end(); ++j) {
    aux_atoms_.push_back(shared_ptr<Atom>(new Atom(**j, displ)));
  }

  common_init1();
  common_init2(false, overlap_thresh_);
}

/************************************************************
*  Merge info from multiple geometries to make one          *
*  supergeometry                                            *
************************************************************/
Geometry::Geometry(vector<shared_ptr<const Geometry> > nmer) :
   spherical_(nmer.front()->spherical_), symmetry_(nmer.front()->symmetry_), schwarz_thresh_(nmer.front()->schwarz_thresh_),
   overlap_thresh_(nmer.front()->overlap_thresh_), external_(nmer.front()->external_)
{
   /************************************************************
   * Going down the list of protected variables, merge the     *
   * data, pick the best ones, or make sure they all match     *
   ************************************************************/
   /* spherical_ should match across the vector*/
   for(auto inmer = nmer.begin(); inmer != nmer.end(); ++inmer) {
      if (spherical_ != (*inmer)->spherical_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of cartesian and spherical bases");
      }
   }
   /* symmetry_ should match across the vector*/
   for(auto inmer = nmer.begin(); inmer != nmer.end(); ++inmer) {
      if (symmetry_ != (*inmer)->symmetry_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of different symmetries");
      }
   }

   /* external field would hopefully match, but for now, if it doesn't, just disable */
   for(auto inmer = nmer.begin(); inmer != nmer.end(); ++inmer) {
      if(!(equal(external_.begin(), external_.end(), (*inmer)->external_.begin()))){
         external_.clear(); break;
      }
   }

   /* atoms_ and aux_atoms_ can be merged */
   vector<shared_ptr<const Atom> > new_atoms;
   vector<shared_ptr<const Atom> > new_aux_atoms;
   for(auto inmer = nmer.begin(); inmer != nmer.end(); ++inmer) {
      auto iatoms = (*inmer)->atoms();
      auto iaux = (*inmer)->aux_atoms();

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
   for(auto inmer = nmer.begin(); inmer != nmer.end(); ++inmer) {
      schwarz_thresh_ = min(schwarz_thresh_, (*inmer)->schwarz_thresh_);
      overlap_thresh_ = min(overlap_thresh_, (*inmer)->overlap_thresh_);
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

Geometry::Geometry(const vector<RefAtom> atoms, const multimap<string, string> geominfo)
  : spherical_(true), input_(""), atoms_(atoms), lmax_(0) {

  schwarz_thresh_ = read_input<double>(geominfo, "schwarz_thresh", 1.0e-12);
  overlap_thresh_ = read_input<double>(geominfo, "thresh_overlap", 1.0e-8);

  // cartesian or not. Look in the atoms info to find out
  spherical_ = atoms.front()->spherical();
  // basis
  auxfile_ = read_input<string>(geominfo, "df_basis", "");
  // symmetry
  symmetry_ = read_input<string>(geominfo, "symmetry", "c1");

  common_init1();


  print_atoms();

  common_init2(true, overlap_thresh_);

  // static external field
  external_ = vector<double>(3);
  external_[0] = read_input<double>(geominfo, "ex", 0.0);
  external_[1] = read_input<double>(geominfo, "ey", 0.0);
  external_[2] = read_input<double>(geominfo, "ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
}


Geometry::~Geometry() {
}


void Geometry::construct_from_atoms(const vector<shared_ptr<const Atom> > atoms, const multimap<string, string> geominfo){

  schwarz_thresh_ = read_input<double>(geominfo, "schwarz_thresh", 1.0e-12);
  overlap_thresh_ = read_input<double>(geominfo, "thresh_overlap", 1.0e-8);

  // cartesian or not.
  const bool cart = read_input<bool>(geominfo, "cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }

  // basis
  auxfile_ = read_input<string>(geominfo, "df_basis", "");
  // symmetry
  symmetry_ = read_input<string>(geominfo, "symmetry", "c1");

}

double Geometry::compute_nuclear_repulsion() {
  double out = 0.0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    const double c = (*iter)->atom_charge();
    for (auto titer = iter + 1; titer != atoms_.end(); ++titer) {
      const double dist = (*iter)->distance(*titer);
      const double charge = c * (*titer)->atom_charge();
      // nuclear repulsion between dummy atoms are not computed here (as in Molpro)
      if (!(*iter)->dummy() || !(*titer)->dummy())
        out += charge / dist;
    }
  }
  return out;
}


void Geometry::print_atoms() const {
  cout << "  === Geometry ===" << endl << endl;
  cout << "  Symmetry: " << symmetry() << endl;
  cout << endl;

  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) (*iter)->print();
  cout << endl;

}


int Geometry::num_count_ncore_only() const {
  int out = 0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    if ((*iter)->atom_number() >= 2 && (*iter)->atom_number() <= 10) out += 2;
    if ((*iter)->atom_number() > 10) throw logic_error("needs to modify Geometry::count_num_ncore for atoms beyond Ne"); // TODO
  }
  return out;
}


int Geometry::num_count_full_valence_nocc() const {
  int out = 0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    if ((*iter)->atom_number() < 2) out += 1;
    if ((*iter)->atom_number() >= 2 && (*iter)->atom_number() <= 10) out += 5;
    if ((*iter)->atom_number() > 10) throw logic_error("needs to modify Geometry::num_count_full_valence_nocc for atoms beyond Ne"); // TODO
  }
  return out;
};


const vector<double> Geometry::compute_grad_vnuc() const {
  // the derivative of Vnuc
  vector<double> grad(natom()*3);
  fill(grad.begin(), grad.end(), 0.0);
  auto giter = grad.begin();
  for (auto aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter, giter+=3) {
    const double ac = (*aiter)->atom_charge();
    for (auto biter = atoms_.begin(); biter != atoms_.end(); ++biter) {
      if (aiter == biter) continue;
      const array<double,3> displ = (*aiter)->displ(*biter);
      const double c = (*biter)->atom_charge() * ac;
      const double dist = (*aiter)->distance(*biter);
      const double dist3 = dist*dist*dist;
      *(giter+0) += c*displ[0]/dist3;
      *(giter+1) += c*displ[1]/dist3;
      *(giter+2) += c*displ[2]/dist3;
    }
  }
  return grad;
}


void Geometry::merge_obs_aux() {
  aux_merged_ = true;
  atoms_.insert(atoms_.end(), aux_atoms_.begin(), aux_atoms_.end());
  for (auto iter = aux_offsets_.begin(); iter != aux_offsets_.end(); ++iter) {
    for (auto citer = iter->begin(); citer != iter->end(); ++citer) {
      *citer += nbasis_;
    }
  }
  offsets_.insert(offsets_.end(), aux_offsets_.begin(), aux_offsets_.end());
  nbasis_ += naux_;
}


vector<double> Geometry::xyz() const {
  vector<double> out;
  out.reserve(3*natom());
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i) {
    out.push_back((*i)->position(0));
    out.push_back((*i)->position(1));
    out.push_back((*i)->position(2));
  }
  return out;
}

array<double,3> Geometry::charge_center() const {
  array<double,3> out{{0.0, 0.0, 0.0}};
  double sum = 0.0;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i) {
    out[0] += (*i)->atom_charge() * (*i)->position(0);
    out[1] += (*i)->atom_charge() * (*i)->position(1);
    out[2] += (*i)->atom_charge() * (*i)->position(2);
    sum += (*i)->atom_charge();
  }
  out[0] /= sum;
  out[1] /= sum;
  out[2] /= sum;
  return out;
}


vector<double> Geometry::schwarz() const {
  vector<shared_ptr<const Shell> > basis;
  for (auto aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter) {
    const vector<shared_ptr<const Shell> > tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
  }
  const int size = basis.size();

  vector<double> schwarz(size * size);
  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<const Shell> b0 = basis[i0];
    for (int i1 = i0; i1 != size; ++i1) {
      const shared_ptr<const Shell> b1 = basis[i1];

      array<shared_ptr<const Shell>,4> input = {{b1, b0, b1, b0}};
#ifdef LIBINT_INTERFACE
      Libint eribatch(input);
#else
      ERIBatch eribatch(input, 1.0);
#endif
      eribatch.compute();
      const double* eridata = eribatch.data();
      const int datasize = eribatch.data_size();
      double cmax = 0.0;
      for (int xi = 0; xi != datasize; ++xi, ++eridata) {
        const double absed = fabs(*eridata);
        if (absed > cmax) cmax = absed;
      }
      schwarz[i0 * size + i1] = cmax;
      schwarz[i1 * size + i0] = cmax;
    }
  }
  return schwarz;
}


bool Geometry::operator==(const Geometry& o) const {
  bool out = true;
  out &= spherical_ == o.spherical_;
  out &= atoms_.size() == o.atoms_.size();
  out &= aux_atoms_.size() == o.aux_atoms_.size();

  for (auto i = atoms_.begin(), j = o.atoms_.begin(); i != atoms_.end(); ++i, ++j) out &= **i == **j;
  for (auto i = aux_atoms_.begin(), j = o.aux_atoms_.begin(); i != aux_atoms_.end(); ++i, ++j) out &= **i == **j;

  out &= aux_merged_ == o.aux_merged_;
  out &= nbasis_ == o.nbasis_;
  out &= nele_ == o.nele_;
  out &= naux_ == o.naux_;
  out &= lmax_ == o.lmax_;
  out &= aux_lmax_ == o.aux_lmax_;

  return out;
}


// connectivity graph

#include <set>
namespace Bagel {
namespace Geometry {
class Node {
  protected:
    const std::shared_ptr<const Atom> myself_;
    const int num_;
    // in order to avoid cyclic references
    std::list<std::weak_ptr<Node> > connected_;


  public:
    Node(const std::shared_ptr<const Atom> o, const int n) : myself_(o), num_(n) {};
    ~Node() {};

    void add_connected(const std::shared_ptr<Node> i) {
      std::weak_ptr<Node> in = i;
      for (auto iter = connected_.begin(); iter != connected_.end(); ++iter)
        if (iter->lock() == i) throw logic_error("Node::add_connected");
      connected_.push_back(in);
    };

    const std::shared_ptr<const Atom> atom() const { return myself_; };
    int num() const { return num_; };

    std::set<std::shared_ptr<Node> > common_center(const std::shared_ptr<Node> o) const {
      std::set<std::shared_ptr<Node> > out;
      for (auto c = connected_.begin(); c != connected_.end(); ++c) {
        if (c->lock()->connected_with(o)) out.insert(c->lock());
      }
      return out;
    };

    bool connected_with(const std::shared_ptr<Node> o) {
      bool out = false;
      for (auto i = connected_.begin(); i != connected_.end(); ++i) {
        if (i->lock() == o) {
          out = true;
          break;
        }
      }
      return out;
    };

};
} }

using namespace Bagel::Geometry;

array<unique_ptr<double[]>,2> Geometry::compute_internal_coordinate() const {
  cout << "    o Connectivitiy analysis" << endl;

  vector<vector<double> > out;
  const size_t size = natom()*3;

  list<shared_ptr<Node> > nodes;
  int n = 0;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i, ++n) {
    if ((*i)->dummy()) throw runtime_error("haven't thought about internal coordinate with dummy atoms (or gradient in genral)");
    nodes.push_back(shared_ptr<Node>(new Node(*i, n)));
  }

  // first pick up bonds
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    auto j = i;
    for (++j ; j != nodes.end(); ++j) {
      // TODO hardwiring 2.0 is NOT a good practice
      if ((*i)->atom()->distance((*j)->atom()) < 2.0*ang2bohr__) {
        (*i)->add_connected(*j);
        (*j)->add_connected(*i);
        cout << "       bond:  " << setw(6) << (*i)->num() << setw(6) << (*j)->num() << "     bond length" <<
                                    setw(10) << setprecision(4) << (*i)->atom()->distance((*j)->atom()) << " bohr" << endl;
        Quatern<double> ip = (*i)->atom()->position();
        Quatern<double> jp = (*j)->atom()->position();
        jp -= ip;  // jp is a vector from i to j
        jp.normalize();
        vector<double> current(size);
        current[3*(*i)->num()+0] =  jp[1];
        current[3*(*i)->num()+1] =  jp[2];
        current[3*(*i)->num()+2] =  jp[3];
        current[3*(*j)->num()+0] = -jp[1];
        current[3*(*j)->num()+1] = -jp[2];
        current[3*(*j)->num()+2] = -jp[3];
        out.push_back(current);
      }
    }
  }

  // then bond angles A-O-B (A<B)
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    auto j = i;
    for (++j; j != nodes.end(); ++j) {
      std::set<std::shared_ptr<Node> > center = (*i)->common_center(*j);
      for (auto c = center.begin(); c != center.end(); ++c) {
        const double theta = (*c)->atom()->angle((*i)->atom(), (*j)->atom());
        cout << "       angle: " << setw(6) << (*c)->num() << setw(6) << (*i)->num() << setw(6) << (*j)->num() <<
                "     angle" << setw(10) << setprecision(4) << theta << " deg" << endl;
        // I found explicit formulas in http://www.ncsu.edu/chemistry/franzen/public_html/nca/int_coord/int_coord.html (thanking the author)
        // 1=A=i, 2=O=c, 3=B=j
        Quatern<double> op = (*c)->atom()->position();
        Quatern<double> ap = (*i)->atom()->position();
        Quatern<double> bp = (*j)->atom()->position();
        Quatern<double> e21 = ap - op;
        Quatern<double> e23 = bp - op;
        const double r21 = e21.norm();
        const double r23 = e23.norm();
        e21.normalize();
        e23.normalize();
        const double rad = theta/rad2deg__;
        Quatern<double> st1 = (e21 * ::cos(rad) - e23) / (r21 * ::sin(rad));
        Quatern<double> st3 = (e23 * ::cos(rad) - e21) / (r23 * ::sin(rad));
        Quatern<double> st2 = (st1 + st3) * (-1.0);
        vector<double> current(size);
        for (int ic = 0; ic != 3; ++ic) {
          current[3*(*i)->num() + ic] = st1[ic+1];
          current[3*(*j)->num() + ic] = st3[ic+1];
          current[3*(*c)->num() + ic] = st2[ic+1];
        }
        out.push_back(current);
      }
    }
  }

  // then dihedral angle (i-c-j-k)
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    for (auto j = nodes.begin(); j != nodes.end(); ++j) {
      if (*i == *j) continue;
      std::set<std::shared_ptr<Node> > center = (*i)->common_center(*j);
      for (auto k = nodes.begin(); k != nodes.end(); ++k) {
        if (!(*k)->connected_with(*j)) continue;
        for (auto c = center.begin(); c != center.end(); ++c) {
          if (*c == *k || *k == *i) continue;
          cout << "    dihedral: " << setw(6) << (*i)->num() << setw(6) << (*c)->num() << setw(6) << (*j)->num() << setw(6) << (*k)->num() <<
                  "     angle" << setw(10) << setprecision(4) << (*c)->atom()->dihedral_angle((*i)->atom(), (*j)->atom(), (*k)->atom()) << " deg" << endl;
          // following J. Molec. Spec. 44, 599 (1972)
          // a=i, b=c, c=j, d=k
          Quatern<double> ap = (*i)->atom()->position();
          Quatern<double> bp = (*c)->atom()->position();
          Quatern<double> cp = (*j)->atom()->position();
          Quatern<double> dp = (*k)->atom()->position();
          Quatern<double> eab = bp - ap;
          Quatern<double> ebc = cp - bp;
          Quatern<double> edc = cp - dp;
          Quatern<double> ecb = bp - cp;
          const double rab = eab.norm();
          const double rbc = ebc.norm();
          const double rcd = edc.norm();
          eab.normalize();
          ebc.normalize();
          edc.normalize();
          ecb.normalize();
          Quatern<double> rotabc = (eab*(-1.0)) * ebc; rotabc[0] = 0.0;
          Quatern<double> rotbcd = (edc*(-1.0)) * ecb; rotbcd[0] = 0.0;
          const double tabc = ::atan2(rotabc.norm(), -eab.ddot(ebc));
          const double tbcd = ::atan2(rotbcd.norm(), -edc.ddot(ecb));

          Quatern<double> sa = (eab * ebc) / (-rab*::pow(::sin(tabc), 2.0));
          Quatern<double> sd = (edc * ecb) / (-rcd*::pow(::sin(tbcd), 2.0));
          Quatern<double> sb = (eab * ebc) * ((rbc-rab*::cos(tabc)) / (rab*rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * (::cos(tbcd) / (rbc*::pow(::sin(tbcd), 2.0)));
          Quatern<double> sc = (eab * ebc) * (::cos(tabc) / (rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * ((rbc-rcd*::cos(tbcd)) / (rcd*rbc*::pow(::sin(tbcd), 2.0)));
          vector<double> current(size);
          for (int ic = 0; ic != 3; ++ic) {
            current[3*(*i)->num() + ic] = sa[ic+1];
            current[3*(*c)->num() + ic] = sb[ic+1];
            current[3*(*j)->num() + ic] = sc[ic+1];
            current[3*(*k)->num() + ic] = sd[ic+1];
          }
          out.push_back(current);
        }
      }
    }
  }

  // debug output
  const int primsize = out.size();
  const int cartsize = 3*natom();

  unique_ptr<double[]> bdag(new double[primsize*cartsize]);
  double* biter = bdag.get();
  for (auto i = out.begin(); i != out.end(); ++i, biter += cartsize)
    copy(i->begin(), i->end(), biter);

  unique_ptr<double[]> bb(new double[primsize*primsize]);
  unique_ptr<double[]> eig(new double[primsize]);
  dgemm_("T", "N", primsize, primsize, cartsize, -1.0, bdag.get(), cartsize, bdag.get(), cartsize, 0.0, bb.get(), primsize);
  const int lwork = primsize*5;
  unique_ptr<double[]> work(new double[lwork]);
  int info;
  dsyev_("V", "U", primsize, bb.get(), primsize, eig.get(), work.get(), lwork, info);
  if (info) throw runtime_error("DSYEV failed in Geometry::compute_internal_coordinate");

  int ninternal = 0;
  for (int i = 0; i != primsize; ++i) {
    if (fabs(eig[i]) > numerical_zero__*::pow(10.0,5)) {
      ++ninternal;
      eig[i] *= -1.0;
    } else {
      break;
    }
  }
  cout << "      Nonredundant internal coordinate generated (dim = " << ninternal << ")" << endl;
  if (ninternal != cartsize-6)
    cout << "       ** caution **  the dimention of internal coordinates is not the same as 3*natom" << endl;

  // form B = U^+ Bprim
  unique_ptr<double[]> bnew(new double[cartsize*cartsize]);
  fill(bnew.get(), bnew.get()+cartsize*cartsize, 0.0);
  dgemm_("T", "T", ninternal, cartsize, primsize, 1.0, bb, primsize, bdag, cartsize, 0.0, bnew, cartsize);

  // form (B^+)^-1 = (BB^+)^-1 B = Lambda^-1 B
  unique_ptr<double[]> bdmnew(new double[cartsize*cartsize]);
  fill(bdmnew.get(), bdmnew.get()+cartsize*cartsize, 0.0);
  for (int j = 0; j != cartsize; ++j)
    for (int i = 0; i != ninternal; ++i)
      bdmnew[i+j*cartsize] = bnew[i+j*cartsize] / eig[i];

  // for debug
  if (false) {
    unique_ptr<double[]> tmp(new double[cartsize*cartsize]);
    fill(tmp.get(), tmp.get()+cartsize*cartsize, 0.0);
    dgemm_("N", "T", cartsize, cartsize, cartsize, 1.0, bdmnew, cartsize, bnew, cartsize, 0.0, tmp, cartsize);
    for (int i = 0; i != cartsize; ++i) {
      for (int j = 0; j != cartsize; ++j) {
        cout << setw(8) << setprecision(3) << tmp[j+i*cartsize];
      }
      cout << endl;
    }
    cout << "====" << endl;
    dgemm_("N", "T", cartsize, cartsize, cartsize, 1.0, bnew, cartsize, bdmnew, cartsize, 0.0, tmp, cartsize);
    for (int i = 0; i != cartsize; ++i) {
      for (int j = 0; j != cartsize; ++j) {
        cout << setw(8) << setprecision(3) << tmp[j+i*cartsize];
      }
      cout << endl;
    }
  }

  return array<unique_ptr<double[]>,2>{{move(bnew), move(bdmnew)}};
}
