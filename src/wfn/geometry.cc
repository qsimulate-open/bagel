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
#include <src/rysint/eribatch.h>
#include <src/rysint/smalleribatch.h>
#include <src/rysint/mixederibatch.h>
#include <src/wfn/geometry.h>
#include <src/io/moldenin.h>
#include <src/util/atommap.h>
#include <src/util/constants.h>
#include <src/util/quatern.h>
#include <src/rysint/libint.h>
#include <src/util/lexical_cast.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap_;

Geometry::Geometry(const boost::property_tree::ptree& geominfo)
  : spherical_(true), input_(""), lmax_(0) {

  schwarz_thresh_ = geominfo.get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo.get<double>("thresh_overlap", 1.0e-8);

  // symmetry
  symmetry_ = geominfo.get<string>("symmetry", "c1");
  transform(symmetry_.begin(), symmetry_.end(), symmetry_.begin(), ::tolower);

  // cartesian or not.
  const bool cart = geominfo.get<bool>("cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }

  /* Set up atoms_ */
  basisfile_ = geominfo.get<string>("basis", "");
  transform(basisfile_.begin(), basisfile_.end(), basisfile_.begin(), ::tolower);
  if (basisfile_ == "") throw runtime_error("There is no basis specification");
  else if (basisfile_ == "molden") {
    string molden_file = geominfo.get<string>("molden_file", "");
    if (molden_file == "") throw runtime_error("No molden_in file provided");

    MoldenIn mfs(molden_file, spherical_);
    mfs.read();
    mfs >> atoms_;
    mfs.close();
  }
  else {
    const bool angstrom = geominfo.get<bool>("angstrom", false);

    auto atoms = geominfo.get_child("geometry");
    for (auto iter = atoms.begin(); iter != atoms.end(); ++iter) {
      string aname = iter->second.get<string>("atom");
      transform(aname.begin(), aname.end(), aname.begin(), ::tolower);
      auto xyz = iter->second.get_child("xyz");
      if (xyz.size() < 3)
        throw runtime_error("XYZ coordinate was not specified correctly");
      array<double,3> positions;
      auto p = positions.begin();
      const double prefac = angstrom ? ang2bohr__ : 1.0;
      for (auto& i : xyz)
        *p++ = lexical_cast<double>(i.second.data()) * prefac;

      if (aname != "q") {
        atoms_.push_back(make_shared<const Atom>(spherical_, aname, positions, basisfile_));
      } else {
        if (symmetry_ != "c1")
          throw runtime_error("External point charges are only allowed in C1 calculations so far.");
        const double charge = lexical_cast<double>(iter->second.get<double>("charge"));
        atoms_.push_back(make_shared<const Atom>(spherical_, aname, positions, charge));
      }

    }
  }
  if(atoms_.empty()) throw runtime_error("No atoms specified at all");

  /* Set up aux_atoms_ */
  auxfile_ = geominfo.get<string>("df_basis", "");
  transform(auxfile_.begin(), auxfile_.end(), auxfile_.begin(), ::tolower);
  if (!auxfile_.empty()) {
    for(auto& iatom : atoms_) {
      if (!iatom->dummy()) {
        aux_atoms_.push_back(make_shared<const Atom>(spherical_, iatom->name(), iatom->position(), auxfile_));
      } else {
        // we need a dummy atom here to be consistent in gradient computations
        aux_atoms_.push_back(iatom);
      }
    }
  }

  // Misc
  aux_merged_ = false;

  common_init1();

  print_atoms();

  common_init2(true, overlap_thresh_);

  // static external field
  external_[0] = geominfo.get<double>("ex", 0.0);
  external_[1] = geominfo.get<double>("ey", 0.0);
  external_[2] = geominfo.get<double>("ez", 0.0);
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

  for (auto& catom : atoms_) {
    nele_ += atommap_.atom_number(catom->name());

    int cc = 0;
    vector<int> coffsets;
    for (auto& it : catom->shells()) {
      coffsets.push_back(nbasis_ + cc);
      const int ang = it->angular_number();
      const int angsize = spherical_ ? (2*ang+1) : (ang+1)*(ang+2)/2;
      cc += angsize * it->num_contracted();
    }
    lmax_ = max(lmax_, catom->lmax());
    nbasis_ += catom->nbasis();
    offsets_.push_back(coffsets);
  }

  if (!auxfile_.empty()) {
    for (auto& catom : aux_atoms_) {
      int cc = 0;
      vector<int> coffsets;
      for (auto& it : catom->shells()) {
        coffsets.push_back(naux_ + cc);
        const int ang = it->angular_number();
        const int angsize = spherical_ ? (2*ang+1) : (ang+1)*(ang+2)/2;
        cc += angsize * it->num_contracted();
      }
      aux_lmax_ = max(lmax_, catom->lmax());
      naux_ += catom->nbasis();
      aux_offsets_.push_back(coffsets);
    }
  }
}


void Geometry::common_init2(const bool print, const double thresh, const bool nodf) {
  // symmetry set-up
  plist_ = make_shared<Petite>(atoms_, symmetry_);
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
    Timer timer;
#ifdef LIBINT_INTERFACE
    df_ = form_fit<DFDist_ints<Libint>>(thresh, true); // true means we construct J^-1/2
#else
    df_ = form_fit<DFDist_ints<ERIBatch>>(thresh, true); // true means we construct J^-1/2
#endif
    cout << "        elapsed time:  " << setw(10) << setprecision(2) << timer.tick() << " sec." << endl << endl;
  }

}



// suitable for geometry updates in optimization
Geometry::Geometry(const Geometry& o, const shared_ptr<const Matrix> displ, const boost::property_tree::ptree& geominfo, const bool rotate, const bool nodf)
  : spherical_(o.spherical_), input_(o.input_), aux_merged_(o.aux_merged_), basisfile_(o.basisfile_),
    auxfile_(o.auxfile_), symmetry_(o.symmetry_), schwarz_thresh_(o.schwarz_thresh_), external_(o.external_), gamma_(o.gamma_) {

  // first construct atoms using displacements
  int iat = 0;
  for (auto i = o.atoms_.begin(), j = o.aux_atoms_.begin(); i != o.atoms_.end(); ++i, ++j, ++iat) {
    array<double,3> cdispl = {{displ->element(0,iat), displ->element(1,iat), displ->element(2,iat)}};
    atoms_.push_back(make_shared<Atom>(**i, cdispl));
    aux_atoms_.push_back(make_shared<Atom>(**j, cdispl));
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
    vector<shared_ptr<const Atom>> newatoms;
    vector<shared_ptr<const Atom>> newauxatoms;
    for (auto i = atoms_.begin(), j = aux_atoms_.begin(); i != atoms_.end(); ++i, ++j) {
      assert((*i)->position() == (*j)->position());
      Quatern<double> source = (*i)->position();
      Quatern<double> target = op * (source - mc) * opd + oc;
      array<double,3> cdispl = (target - source).ijk();

      newatoms.push_back(make_shared<Atom>(**i, cdispl));
      newauxatoms.push_back(make_shared<Atom>(**j, cdispl));
    }
    atoms_ = newatoms;
    aux_atoms_ = newauxatoms;
  }

  common_init1();
  overlap_thresh_ = geominfo.get<double>("thresh_overlap", 1.0e-8);
  common_init2(false, overlap_thresh_, nodf);
}


Geometry::Geometry(const Geometry& o, const array<double,3> displ)
  : spherical_(o.spherical_), input_(o.input_), aux_merged_(o.aux_merged_), basisfile_(o.basisfile_),
    auxfile_(o.auxfile_), symmetry_(o.symmetry_), schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_),
    external_(o.external_), gamma_(o.gamma_) {

  // first construct atoms using displacements

  for (auto& i : o.atoms_) {
    atoms_.push_back(make_shared<Atom>(*i, displ));
  }

  for (auto& j : o.aux_atoms_) {
    aux_atoms_.push_back(make_shared<Atom>(*j, displ));
  }

  common_init1();
  common_init2(false, overlap_thresh_);
}

/************************************************************
*  Merge info from multiple geometries to make one          *
*  supergeometry                                            *
************************************************************/
Geometry::Geometry(vector<shared_ptr<const Geometry>> nmer) :
   spherical_(nmer.front()->spherical_), symmetry_(nmer.front()->symmetry_), schwarz_thresh_(nmer.front()->schwarz_thresh_),
   overlap_thresh_(nmer.front()->overlap_thresh_), external_(nmer.front()->external_)
{
   /************************************************************
   * Going down the list of protected variables, merge the     *
   * data, pick the best ones, or make sure they all match     *
   ************************************************************/
   /* spherical_ should match across the vector*/
   for(auto& inmer : nmer) {
      if (spherical_ != (inmer)->spherical_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of cartesian and spherical bases");
      }
   }

   /* symmetry_ should match across the vector*/
   for(auto& inmer : nmer) {
      if (symmetry_ != (inmer)->symmetry_) {
         throw runtime_error("Attempting to construct a geometry that is a mixture of different symmetries");
      }
   }

   /* external field would hopefully match, but for now, if it doesn't, just disable */
   for(auto& inmer : nmer) {
      if(!(equal(external_.begin(), external_.end(), inmer->external_.begin()))){
         fill(external_.begin(), external_.end(), 0.0); break;
      }
   }

   /* atoms_ and aux_atoms_ can be merged */
   vector<shared_ptr<const Atom>> new_atoms;
   vector<shared_ptr<const Atom>> new_aux_atoms;
   for(auto& inmer : nmer) {
      auto iatoms = inmer->atoms();
      auto iaux = inmer->aux_atoms();

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
   for(auto& inmer : nmer) {
      schwarz_thresh_ = min(schwarz_thresh_, inmer->schwarz_thresh_);
      overlap_thresh_ = min(overlap_thresh_, inmer->overlap_thresh_);
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

Geometry::Geometry(const vector<shared_ptr<const Atom>> atoms, const boost::property_tree::ptree& geominfo)
  : spherical_(true), input_(""), atoms_(atoms), lmax_(0) {

  schwarz_thresh_ = geominfo.get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo.get<double>("thresh_overlap", 1.0e-8);

  // cartesian or not. Look in the atoms info to find out
  spherical_ = atoms.front()->spherical();
  // basis
  auxfile_ = geominfo.get<string>("df_basis", "");
  if (!auxfile_.empty()) {
    if (!aux_atoms_.empty()) throw logic_error("programming error in the Geometry constructor with vector<shared_ptr<Atom>>");
    for (auto& i : atoms_)
      aux_atoms_.push_back(make_shared<const Atom>(i->spherical(), i->name(), i->position(), auxfile_));
  }
  // symmetry
  symmetry_ = geominfo.get<string>("symmetry", "c1");

  common_init1();

  print_atoms();

  common_init2(true, overlap_thresh_);

  // static external field
  external_[0] = geominfo.get<double>("ex", 0.0);
  external_[1] = geominfo.get<double>("ey", 0.0);
  external_[2] = geominfo.get<double>("ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
}


void Geometry::construct_from_atoms(const vector<shared_ptr<const Atom>> atoms, const boost::property_tree::ptree& geominfo){

  schwarz_thresh_ = geominfo.get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo.get<double>("thresh_overlap", 1.0e-8);

  // cartesian or not.
  const bool cart = geominfo.get<bool>("cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }

  // basis
  auxfile_ = geominfo.get<string>("df_basis", "");
  // symmetry
  symmetry_ = geominfo.get<string>("symmetry", "c1");

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
  cout << "  *** Geometry ***" << endl << endl;
  cout << "  Symmetry: " << symmetry() << endl;
  cout << endl;

  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) (*iter)->print();
  cout << endl;

}


int Geometry::num_count_ncore_only() const {
  int out = 0;
  for (auto& it : atoms_) {
    if (it->atom_number() >= 2 && it->atom_number() <= 10) out += 2;
    if (it->atom_number() > 10) throw logic_error("needs to modify Geometry::count_num_ncore for atoms beyond Ne"); // TODO
  }
  return out;
}


int Geometry::num_count_full_valence_nocc() const {
  int out = 0;
  for (auto& it : atoms_) {
    if (it->atom_number() < 2) out += 1;
    if (it->atom_number() >= 2 && it->atom_number() <= 10) out += 5;
    if (it->atom_number() > 10) throw logic_error("needs to modify Geometry::num_count_full_valence_nocc for atoms beyond Ne"); // TODO
  }
  return out;
};


const shared_ptr<const Matrix> Geometry::compute_grad_vnuc() const {
  // the derivative of Vnuc
  auto grad = make_shared<Matrix>(3, natom());
  int i = 0;
  for (auto& a : atoms_) {
    const double ac = a->atom_charge();
    for (auto& b : atoms_) {
      if (a == b) continue;
      const array<double,3> displ = a->displ(b);
      const double c = b->atom_charge() * ac;
      const double dist = a->distance(b);
      const double dist3 = dist*dist*dist;
      grad->element(0,i) += c*displ[0]/dist3;
      grad->element(1,i) += c*displ[1]/dist3;
      grad->element(2,i) += c*displ[2]/dist3;
    }
    ++i;
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


shared_ptr<const Matrix> Geometry::xyz() const {
  auto out = make_shared<Matrix>(3, natom());
  int iat = 0;
  for (auto& i : atoms_) {
    out->element(0, iat) = i->position(0);
    out->element(1, iat) = i->position(1);
    out->element(2, iat) = i->position(2);
    ++iat;
  }
  return out;
}


array<double,3> Geometry::charge_center() const {
  array<double,3> out{{0.0, 0.0, 0.0}};
  double sum = 0.0;
  for (auto& i : atoms_) {
    out[0] += i->atom_charge() * i->position(0);
    out[1] += i->atom_charge() * i->position(1);
    out[2] += i->atom_charge() * i->position(2);
    sum += i->atom_charge();
  }
  out[0] /= sum;
  out[1] /= sum;
  out[2] /= sum;
  return out;
}


vector<double> Geometry::schwarz() const {
  vector<shared_ptr<const Shell>> basis;
  for (auto aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter) {
    const vector<shared_ptr<const Shell>> tmp = (*aiter)->shells();
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
namespace bagel {
namespace geometry {
class Node {
  protected:
    const std::shared_ptr<const Atom> myself_;
    const int num_;
    // in order to avoid cyclic references
    std::list<std::weak_ptr<Node>> connected_;


  public:
    Node(const std::shared_ptr<const Atom> o, const int n) : myself_(o), num_(n) {};
    ~Node() {};

    void add_connected(const std::shared_ptr<Node> i) {
      std::weak_ptr<Node> in = i;
      for (auto& iter : connected_)
        if (iter.lock() == i) throw logic_error("Node::add_connected");
      connected_.push_back(in);
    };

    const std::shared_ptr<const Atom> atom() const { return myself_; };
    int num() const { return num_; };

    std::set<std::shared_ptr<Node>> common_center(const std::shared_ptr<Node> o) const {
      std::set<std::shared_ptr<Node>> out;
      for (auto& c : connected_) {
        if (c.lock()->connected_with(o)) out.insert(c.lock());
      }
      return out;
    };

    bool connected_with(const std::shared_ptr<Node> o) {
      bool out = false;
      for (auto& i : connected_) {
        if (i.lock() == o) {
          out = true;
          break;
        }
      }
      return out;
    };

};
} }

using namespace bagel::geometry;

array<unique_ptr<double[]>,2> Geometry::compute_internal_coordinate() const {
  cout << "    o Connectivitiy analysis" << endl;

  vector<vector<double>> out;
  const size_t size = natom()*3;

  list<shared_ptr<Node>> nodes;
  int n = 0;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i, ++n) {
    if ((*i)->dummy()) throw runtime_error("haven't thought about internal coordinate with dummy atoms (or gradient in genral)");
    nodes.push_back(make_shared<Node>(*i, n));
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
      std::set<std::shared_ptr<Node>> center = (*i)->common_center(*j);
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
      std::set<std::shared_ptr<Node>> center = (*i)->common_center(*j);
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

  auto bb = make_shared<Matrix>(primsize,primsize,true);
  dgemm_("T", "N", primsize, primsize, cartsize, -1.0, bdag.get(), cartsize, bdag.get(), cartsize, 0.0, bb->data(), primsize);
  unique_ptr<double[]> eig(new double[primsize]);
  bb->diagonalize(eig.get());

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
  dgemm_("T", "T", ninternal, cartsize, primsize, 1.0, bb->data(), primsize, bdag.get(), cartsize, 0.0, bnew.get(), cartsize);

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


shared_ptr<const Geometry> Geometry::relativistic(const bool do_gaunt) const {
  cout << "  *** Geometry (Relativistic) ***" << endl;
  Timer timer;
  // basically the same
  auto geom = make_shared<Geometry>(*this);

  // except for atoms_->shells
  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->relativistic());
  geom->atoms_ = atom;
  geom->df_->average_3index();
  geom->dfs_  = geom->form_fit<DFDist_ints<SmallERIBatch>>(overlap_thresh_, true, 0.0, true);
  if (do_gaunt)
    geom->dfsl_ = geom->form_fit<DFDist_ints<MixedERIBatch>>(overlap_thresh_, true, 0.0, true);

  // suppress some of the printing
  resources__->proc()->set_print_level(2);
  cout << endl;
  timer.tick_print("Geometry relativistic (total)");
  cout << endl;
  return geom;
}


void Geometry::discard_relativistic() const {
  dfs_ = shared_ptr<DFDist>();
  dfsl_ = shared_ptr<DFDist>();
}
