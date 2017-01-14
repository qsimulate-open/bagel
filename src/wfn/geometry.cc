//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: geometry.cc
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


#include <src/wfn/geometry.h>
#include <src/df/complexdf.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/smalleribatch.h>
#include <src/integral/rys/mixederibatch.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/comprys/complexsmalleribatch.h>
#include <src/integral/comprys/complexmixederibatch.h>
#include <src/integral/libint/libint.h>
#include <src/util/io/moldenin.h>
#include <src/util/math/quatern.h>

using namespace std;
using namespace bagel;


BOOST_CLASS_EXPORT_IMPLEMENT(Geometry)

Geometry::Geometry(shared_ptr<const PTree> geominfo) : magnetism_(false), do_periodic_df_(false) {

  // members of Molecule
  spherical_ = true;
  lmax_ = 0;
  dofmm_   = geominfo->get<bool>("cfmm", false);

  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);

  // symmetry
  symmetry_ = to_lower(geominfo->get<string>("symmetry", "c1"));

  // cartesian or not.
  const bool cart = geominfo->get<bool>("cartesian", false);
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }
  const bool angstrom = geominfo->get<bool>("angstrom", false);

  /* Set up lattice */
  shared_ptr<const PTree> vectors = geominfo->get_child_optional("primitive_vectors");
  if (vectors) {
    int dim = 0;
    for (auto& ivec : *vectors) {
      string id = "a" + to_string(dim+1);
      array<double, 3> vec = ivec->get_array<double, 3>(id);
      if (angstrom) for (auto& i : vec) i /= angstrom ? au2angstrom__ : 1.0;
      primitive_vectors_.push_back(vec);
      ++dim;
    }
  }

  // DKH2 Hamiltonian
  dkh_ = geominfo->get<bool>("dkh", false);

  // static external magnetic field
  magnetic_field_ = geominfo->get_array<double,3>("magnetic_field", {{0.0, 0.0, 0.0}});
  const bool tesla = geominfo->get<bool>("tesla", false);
  if (tesla)
    for (int i = 0; i != 3; ++i)
      magnetic_field_[i] /= au2tesla__;
  set_london(geominfo);

  /* Set up atoms_ */
  basisfile_ = geominfo->get<string>("basis", "");
  use_finite_ = geominfo->get<bool>("finite_nucleus", false);
  if (basisfile_ == "") {
    throw runtime_error("There is no basis specification");
  } else if (basisfile_ == "molden") {
    string molden_file = geominfo->get<string>("molden_file", "");
    if (molden_file == "") throw runtime_error("No molden_in file provided");

    MoldenIn mfs(molden_file, spherical_);
    mfs.read();
    mfs >> atoms_;
  } else {

    // read the default basis file
    shared_ptr<const PTree> bdata = PTree::read_basis(basisfile_);
    shared_ptr<const PTree> elem = geominfo->get_child_optional("_basis");

    auto atoms = geominfo->get_child("geometry");
    use_ecp_basis_ = (basisfile_.find("ecp") != string::npos) ? true : false;
    for (auto& a : *atoms)
      atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem, false, use_ecp_basis_, use_finite_));
  }
  if (atoms_.empty()) throw runtime_error("No atoms specified at all");
  for (auto& i : atoms_)
    if (symmetry_ != "c1" && i->dummy())
      throw runtime_error("External point charges are only allowed in C1 calculations so far.");

  /* Set up aux_atoms_ */
  auxfile_ = geominfo->get<string>("df_basis", "");  // default value for non-DF HF.
  if (!auxfile_.empty()) {
    if (!primitive_vectors_.empty()) do_periodic_df_ = true;
    // read the default aux basis file
    shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (basisfile_ == "molden") {
      for(auto& iatom : atoms_) {
        if (!iatom->dummy()) {
          aux_atoms_.push_back(make_shared<const Atom>(spherical_, iatom->name(), iatom->position(), auxfile_, make_pair(auxfile_, bdata), elem));
        } else {
          // we need a dummy atom here to be consistent in gradient computations
          aux_atoms_.push_back(iatom);
        }
      }
    } else {
      auto atoms = geominfo->get_child("geometry");
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem, true));
    }
  }

  // Misc
  aux_merged_ = false;

  common_init1();

  print_atoms();

  common_init2(true, overlap_thresh_);
  get_electric_field(geominfo);
}


void Geometry::common_init2(const bool print, const double thresh, const bool nodf) {

  if (london_ || nonzero_magnetic_field()) init_magnetism();

  if (!auxfile_.empty() && !nodf && !do_periodic_df_ && !dofmm_) {
    if (print) cout << "  Number of auxiliary basis functions: " << setw(8) << naux() << endl << endl;
    cout << "  Since a DF basis is specified, we compute 2- and 3-index integrals:" << endl;
    const double scale = magnetism_ ? 2.0 : 1.0;
    cout << "    o Being stored without compression. Storage requirement is "
         << setprecision(3) << static_cast<size_t>(naux_)*nbasis()*nbasis()*scale*8.e-9 << " GB" << endl;
    Timer timer;
    compute_integrals(thresh);
    cout << "        elapsed time:  " << setw(10) << setprecision(2) << timer.tick() << " sec." << endl << endl;
  }

  // symmetry set-up
  plist_ = make_shared<Petite>(atoms_, symmetry_);
  nirrep_ = plist_->nirrep();

  if (print) {
    cout << endl;
    cout << "  Number of basis functions: " << setw(8) << nbasis() << endl;
    cout << "  Number of electrons      : " << setw(8) << nele() << endl << endl;
  }

  nuclear_repulsion_ = compute_nuclear_repulsion();
  get_shellpairs();

  assert(magnetism_ ? (london_ || nonzero_magnetic_field()) : (!london_ && !nonzero_magnetic_field()));
}


// suitable for geometry updates in optimization
Geometry::Geometry(const Geometry& o, shared_ptr<const Matrix> displ, shared_ptr<const PTree> geominfo, const bool rotate, const bool nodf)
  : schwarz_thresh_(o.schwarz_thresh_), dkh_(o.dkh_), magnetism_(false), london_(o.london_), use_finite_(o.use_finite_), use_ecp_basis_(o.use_ecp_basis_),
    do_periodic_df_(o.do_periodic_df_) {

  // Members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  gamma_ = o.gamma_;
  external_ = o.external_;
  magnetic_field_ = o.magnetic_field_;
  dofmm_ = o.dofmm_;

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
    int iatom = 0;
    for ( ; iatom != natom(); ++iatom) {
      Quatern<double> oa = o.atoms(iatom)->position();
      Quatern<double> ma =   atoms(iatom)->position();
      Quatern<double> od = oa - oc;
      Quatern<double> md = ma - mc;
      // if the charge center coincide with the location of the atom, skip
      if (od.norm() < 0.1 || md.norm() < 0.1)
        continue;
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
      break;
    }

    // (3) plane of center of charges, first and second atoms.
    if (natom() > 2) {
      assert(natom() == o.natom());
      for (int jatom = 0; jatom != natom(); ++jatom) {
        if (iatom == jatom) continue;
        Quatern<double> oa0 = o.atoms(iatom)->position();
        Quatern<double> ma0 =   atoms(iatom)->position();
        Quatern<double> oa1 = o.atoms(jatom)->position();
        Quatern<double> ma1 =   atoms(jatom)->position();
        Quatern<double> mc = charge_center();
        Quatern<double> od = (oa0 - oc) * (oa1 - oc);
        Quatern<double> md = (ma0 - mc) * (ma1 - mc);
        od[0] = 0.0;
        md[0] = 0.0;
        if (od.norm() < 1.0e-5 || md.norm() < 1.0e-5)
          continue;
        od.normalize();
        md.normalize();
        Quatern<double> op = md * od;
        op[0] = 1.0 - op[0];
        op.normalize();
        Quatern<double> opd = op.dagger();

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
        break;
      }
    }
  }

  common_init1();
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);
  set_london(geominfo);
  common_init2(false, overlap_thresh_, nodf);
  if (o.magnetism())
    throw logic_error("Geometry optimization in a magnetic field has not been set up or verified; use caution.");
}

Geometry::Geometry(const Geometry& o, const array<double,3> displ)
  : schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_), dkh_(o.dkh_), magnetism_(false),
    london_(o.london_), use_finite_(o.use_finite_), use_ecp_basis_(o.use_ecp_basis_), do_periodic_df_(o.do_periodic_df_) {

  // members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  gamma_ = o.gamma_;
  external_ = o.external_;
  magnetic_field_ = o.magnetic_field_;
  dofmm_ = o.dofmm_;

  // first construct atoms using displacements
  for (auto& i : o.atoms_) {
    atoms_.push_back(make_shared<Atom>(*i, displ));
  }

  for (auto& j : o.aux_atoms_) {
    aux_atoms_.push_back(make_shared<Atom>(*j, displ));
  }

  common_init1();
  common_init2(false, overlap_thresh_);
  if (o.magnetism())
    throw logic_error("Geometry displacement in a magnetic field has not been set up or verified; use caution.");
}


// used when a new Geometry block is provided in input
Geometry::Geometry(const Geometry& o, shared_ptr<const PTree> geominfo, const bool discard)
  : schwarz_thresh_(o.schwarz_thresh_), overlap_thresh_(o.overlap_thresh_), magnetism_(false),
    london_(o.london_), use_finite_(o.use_finite_), use_ecp_basis_(o.use_ecp_basis_), do_periodic_df_(o.do_periodic_df_) {

  // members of Molecule
  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  symmetry_ = o.symmetry_;
  external_ = o.external_;
  atoms_ = o.atoms_;
  aux_atoms_ = o.aux_atoms_;
  gamma_ = o.gamma_;
  magnetic_field_ = o.magnetic_field_;
  dofmm_ = o.dofmm_;

  // check all the options
  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", schwarz_thresh_);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", overlap_thresh_);
  symmetry_ = to_lower(geominfo->get<string>("symmetry", symmetry_));

  spherical_ = !geominfo->get<bool>("cartesian", !spherical_);
  dofmm_ = geominfo->get<bool>("cfmm", false);
  dkh_ = geominfo->get<bool>("dkh", dkh_);

  // check if a magnetic field has been supplied
  auto newfield = geominfo->get_child_optional("magnetic_field");
  if (newfield) {
    magnetic_field_ = geominfo->get_array<double,3>("magnetic_field");
    const bool tesla = geominfo->get<bool>("tesla", false);
    if (tesla)
      for (int i = 0; i != 3; ++i)
        magnetic_field_[i] /= au2tesla__;

    const string basis = geominfo->get<string>("basis_type", london_ ? "giao" : "gaussian");
    if (basis == "giao" || basis == "london") london_ = true;
    else if (basis == "gaussian") london_ = false;
    else throw runtime_error("Basis set type not recognized; should be Gaussian or London");
  }

  // check if we need to construct shells and integrals
  auto atoms = geominfo->get_child_optional("geometry");
  const string prevbasis = basisfile_;
  basisfile_ = geominfo->get<string>("basis", basisfile_);
  use_finite_ = geominfo->get<bool>("finite_nucleus", use_finite_);
  // if so, construct atoms
  if (prevbasis != basisfile_ || atoms || newfield) {
    use_ecp_basis_ = (basisfile_.find("ecp") != string::npos) ? true : false;
    atoms_.clear();
    shared_ptr<const PTree> bdata = PTree::read_basis(basisfile_);
    shared_ptr<const PTree> elem = geominfo->get_child_optional("_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(basisfile_, bdata), elem, false, use_ecp_basis_, use_finite_));
    } else {
      for (auto& a : o.atoms_)
        atoms_.push_back(make_shared<const Atom>(*a, spherical_, basisfile_, make_pair(basisfile_, bdata), elem));
    }
  }
  const string prevaux = auxfile_;
  auxfile_ = geominfo->get<string>("df_basis", auxfile_);
  if (prevaux != auxfile_ || atoms) {
    aux_atoms_.clear();
    shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (atoms) {
      const bool angstrom = geominfo->get<bool>("angstrom", false);
      for (auto& a : *atoms)
        aux_atoms_.push_back(make_shared<const Atom>(a, spherical_, angstrom, make_pair(auxfile_, bdata), elem));
    } else {
      for (auto& a : o.atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(*a, spherical_, auxfile_, make_pair(auxfile_, bdata), elem));
    }
  }

  common_init1();

  if (o.basisfile_ != basisfile_ || o.auxfile_ != auxfile_ || atoms || newfield) {
    // discard the previous one before we compute the new one. Note that df_'s are mutable... too bad, I know..
    if (discard)
      o.discard_df();
    common_init2(true, overlap_thresh_);
  } else {
    df_ = o.df_;
    dfs_ = o.dfs_;
    dfsl_ = o.dfsl_;
    common_init2(true, overlap_thresh_, true /* not to calculate integrals */);
  }
  external_ = o.external_;
}


/************************************************************
*  Merge info from multiple geometries to make one          *
*  supergeometry                                            *
************************************************************/
Geometry::Geometry(vector<shared_ptr<const Geometry>> nmer, const bool nodf) :
  schwarz_thresh_(nmer.front()->schwarz_thresh_), overlap_thresh_(nmer.front()->overlap_thresh_), dkh_(nmer.front()->dkh()), magnetism_(false), london_(nmer.front()->london_),
  use_finite_(nmer.front()->use_finite_), use_ecp_basis_(nmer.front()->use_ecp_basis_),  do_periodic_df_(false) {

  // A member of Molecule
  spherical_ = nmer.front()->spherical_;
  symmetry_ = nmer.front()->symmetry_;
  external_ = nmer.front()->external_;
  magnetic_field_ = nmer.front()->magnetic_field_;
  dofmm_ = nmer.front()->dofmm_;

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
  common_init2(true,overlap_thresh_, nodf);

  // static external field
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;

  if (nmer.front()->magnetism())
    throw logic_error("Combining Geometries in a magnetic field has not been set up or verified; use caution.");
}


// used in SCF initial guess.
Geometry::Geometry(const vector<shared_ptr<const Atom>> atoms, shared_ptr<const PTree> geominfo) : magnetism_(false), do_periodic_df_(false) {

  spherical_ = true;
  lmax_ = 0;
  dofmm_ = false;

  atoms_ = atoms;

  auto atomlist = geominfo->get_child_optional("geometry");
  const bool angstrom = geominfo->get<bool>("angstrom", false);

  schwarz_thresh_ = geominfo->get<double>("schwarz_thresh", 1.0e-12);
  overlap_thresh_ = geominfo->get<double>("thresh_overlap", 1.0e-8);

  // cartesian or not. Look in the atoms info to find out
  spherical_ = atoms.front()->spherical();
  // basis
  auxfile_ = geominfo->get<string>("df_basis", "");
  if (!auxfile_.empty()) {
    // read the default basis file
    shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    shared_ptr<const PTree> elem = geominfo->get_child_optional("_df_basis");
    if (atomlist) {
      for (auto& i : *atomlist)
        aux_atoms_.push_back(make_shared<const Atom>(i, spherical_, angstrom, make_pair(auxfile_, bdata), elem, true));
    } else {
      // in the molden case
      for (auto& i : atoms_)
        aux_atoms_.push_back(make_shared<const Atom>(i->spherical(), i->name(), i->position(), auxfile_, make_pair(auxfile_, bdata), elem));
    }
  }
  // symmetry
  symmetry_ = geominfo->get<string>("symmetry", "c1");

  common_init1();

  print_atoms();

  dkh_ = geominfo->get<bool>("dkh", false);

  // static external magnetic field
  magnetic_field_ = geominfo->get_array<double,3>("magnetic_field", {{0.0, 0.0, 0.0}});
  const bool tesla = geominfo->get<bool>("tesla", false);
  if (tesla)
    for (int i=0; i!=3; ++i)
      magnetic_field_[i] /= au2tesla__;
  set_london(geominfo);

  common_init2(true, overlap_thresh_);
  get_electric_field(geominfo);
}


shared_ptr<const Matrix> Geometry::compute_grad_vnuc() const {
  // the derivative of Vnuc
  auto grad = make_shared<Matrix>(3, natom());
  int i = 0;
  for (auto& a : atoms_) {
    const double ac = a->atom_charge();
    if (i % mpi__->size() == mpi__->rank()) {
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
    }
    ++i;
  }
  return grad;
}


void Geometry::get_shellpairs() {

  vector<int> offsets;
  vector<shared_ptr<const Shell>> basis;
  for (int n = 0; n != natom(); ++n) {
    const vector<int> tmpoff = offset(n);
    offsets.insert(offsets.end(), tmpoff.begin(), tmpoff.end());
    const vector<shared_ptr<const Shell>> tmpsh = atoms_[n]->shells();
    basis.insert(basis.end(), tmpsh.begin(), tmpsh.end());
  }
  const int nsh = basis.size();

  shellpairs_.resize(nsh * nsh);
  for (int i0 = 0; i0 != nsh; ++i0) {
    for (int i1 = 0; i1 != nsh; ++i1) {
      const int i01 = i0 * nsh + i1;
      shellpairs_[i01] = make_shared<const ShellPair>(array<shared_ptr<const Shell>, 2>{{basis[i1], basis[i0]}},
                                                      array<int, 2>{{offsets[i1], offsets[i0]}}, make_pair(i1, i0));
    }
  }
}


vector<double> Geometry::schwarz() const {

  const int nsp = shellpairs_.size();
  vector<double> schwarz(nsp);
  for (int i = 0; i != nsp; ++i)
    schwarz[i] = shellpairs_[i]->schwarz();

  return schwarz;
}


void Geometry::get_electric_field(shared_ptr<const PTree>& geominfo) {
  // static external electric field
  external_[0] = geominfo->get<double>("ex", 0.0);
  external_[1] = geominfo->get<double>("ey", 0.0);
  external_[2] = geominfo->get<double>("ez", 0.0);
  if (external())
  cout << "  * applying an external electric field (" << setprecision(3) << setw(7) << external_[0] << ", "
                                                                         << setw(7) << external_[1] << ", "
                                                                         << setw(7) << external_[2] << ") a.u." << endl << endl;
}


shared_ptr<const Geometry> Geometry::periodic(vector<shared_ptr<const Atom>> atoms) const {

  auto out = make_shared<Geometry>(*this);

  vector<shared_ptr<const Atom>> aux_atoms;
  if (!auxfile_.empty()) {
    shared_ptr<const PTree> bdata = PTree::read_basis(auxfile_);
    for (auto& a : atoms)
      aux_atoms.push_back(make_shared<const Atom>(*a, spherical_, auxfile_, make_pair(auxfile_, bdata), nullptr));
  }
  out->atoms_ = atoms;
  out->aux_atoms_ = aux_atoms;
  out->do_periodic_df_ = false;

  out->common_init1();
  out->common_init2(true, overlap_thresh_);

  return make_shared<const Geometry>(*out);
}


shared_ptr<const Geometry> Geometry::relativistic(const bool do_gaunt, const bool do_coulomb) const {
  cout << "  *** Geometry (Relativistic) ***" << endl;
  Timer timer;
  // basically the same
  auto geom = make_shared<Geometry>(*this);

  // except for atoms_->shells
  vector<shared_ptr<const Atom>> atom;

  for (auto& i : atoms_)
    atom.push_back(!magnetism_ ? i->relativistic() : i->relativistic(magnetic_field_, london_));

  geom->atoms_ = atom;

  if (do_coulomb)
    geom->compute_relativistic_integrals(do_gaunt);

  cout << endl;
  timer.tick_print("Geometry relativistic (total)");
  cout << endl;
  return geom;
}


void Geometry::compute_relativistic_integrals(const bool do_gaunt) {
  df_->average_3index();
  shared_ptr<Matrix> d2 = df_->data2()->copy();

  if (!magnetism_) {
    dfs_  = form_fit<DFDist_ints<SmallERIBatch>>(overlap_thresh_, true, 0.0, true, d2);
    if (do_gaunt)
      dfsl_ = form_fit<DFDist_ints<MixedERIBatch>>(overlap_thresh_, true, 0.0, true, d2);
  } else {
    dfs_  = form_fit<ComplexDFDist_ints<ComplexSmallERIBatch>>(overlap_thresh_, true, 0.0, true, d2);
    if (do_gaunt)
      dfsl_ = form_fit<ComplexDFDist_ints<ComplexMixedERIBatch>>(overlap_thresh_, true, 0.0, true, d2);
  }

  // suppress some of the printing
  resources__->proc()->set_print_level(2);
}


void Geometry::discard_relativistic() const {
  dfs_.reset();
  dfsl_.reset();
}


void Geometry::set_london(shared_ptr<const PTree>& geominfo) {
  const string basis_type = to_lower(geominfo->get<string>("basis_type", (nonzero_magnetic_field() ? "london" : "gaussian")));
  if (basis_type == "giao" || basis_type == "london") london_ = true;
  else if (basis_type == "gaussian") london_ = false;
  else
    throw runtime_error("Invalid basis type entered - should be london or gaussian");
}


void Geometry::compute_integrals(const double thresh) const {
#ifdef LIBINT_INTERFACE
  if (!magnetism_)
    df_ = form_fit<DFDist_ints<Libint>>(thresh, true); // true means we construct J^-1/2
#else
  if (!magnetism_)
    df_ = form_fit<DFDist_ints<ERIBatch>>(thresh, true); // true means we construct J^-1/2
#endif
  else
    df_ = form_fit<ComplexDFDist_ints<ComplexERIBatch>>(thresh, true); // true means we construct J^-1/2
}


void Geometry::init_magnetism() {
  magnetism_ = true;

  if (nonzero_magnetic_field()) {
    cout << (london_ ? "  Using London orbital basis to enforce gauge-invariance"
                     : "  Using a common gauge origin - NOT RECOMMENDED for accurate calculations (use a London orbital basis instead).") << endl;

    cout << "  Applied magnetic field:  (" << setprecision(4) << setw(7) << magnetic_field_[0] << ", "
                                                              << setw(7) << magnetic_field_[1] << ", "
                                                              << setw(7) << magnetic_field_[2] << ") a.u." << endl;

    const double fieldsqr = magnetic_field_[0]*magnetic_field_[0] + magnetic_field_[1]*magnetic_field_[1] + magnetic_field_[2]*magnetic_field_[2];
    cout << setprecision(0) << "  Field strength = " << au2tesla__*sqrt(fieldsqr) << " T" << endl << endl;
  } else {
    cout << "  Zero magnetic field - This computation would be more efficient with a standard basis." << endl;
  }

  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->apply_magnetic_field(magnetic_field_, london_));
  atoms_ = atom;
}
