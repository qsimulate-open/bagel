//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molecule.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/molecule/molecule.h>
#include <src/molecule/molecule_connect.h>
#include <src/util/constants.h>
#include <src/util/atommap.h>
#include <src/util/math/quatern.h>

using namespace std;
using namespace bagel;
using molecule_details::Node;
using molecule_details::adf_rho;

const static AtomMap atommap;


BOOST_CLASS_EXPORT_IMPLEMENT(Molecule)


// Initialize Molecule with displacement
Molecule::Molecule(const Molecule& o, shared_ptr<const Matrix> displ, const bool rotate) {

  spherical_ = o.spherical_;
  aux_merged_ = o.aux_merged_;
  basisfile_ = o.basisfile_;
  auxfile_ = o.auxfile_;
  external_ = o.external_;
  magnetic_field_ = o.magnetic_field_;
  skip_self_interaction_ = o.skip_self_interaction_;
  cap_ = o.cap_;

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
}


double Molecule::compute_nuclear_repulsion() {
  double out = 0.0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    double c = (*iter)->atom_charge();
    if ((*iter)->use_ecp_basis()) c -= (*iter)->ecp_parameters()->ecp_ncore();
    for (auto titer = iter + 1; titer != atoms_.end(); ++titer) {
      // nuclear repulsion between dummy atoms are not computed here (as in Molpro)
      if (!skip_self_interaction_ || (!(*iter)->dummy() || !(*titer)->dummy())) {
        const double dist = (*iter)->distance(*titer);
        double tc = (*titer)->atom_charge();
        if ((*titer)->use_ecp_basis()) tc -= (*titer)->ecp_parameters()->ecp_ncore();
        const double charge = c * tc;
// it turned out that the deviation is numerically zero in double precision if the exponent is more than 1.0e-6. I will leave them for a while
#ifdef BEYOND_DOUBLE
        if (!(*iter)->finite_nucleus() && !(*titer)->finite_nucleus()) { // both point charges
#endif
          out += charge / dist;
#ifdef BEYOND_DOUBLE
        } else if ((*iter)->finite_nucleus() && (*titer)->finite_nucleus()) { // both gaussian charges
          const double gamma0 = (*iter)->atom_exponent();
          const double gamma1 = (*titer)->atom_exponent();
          out += charge * erf(sqrt(gamma0 * gamma1 / (gamma0 + gamma1)) * dist) / dist;
        } else { // one point charge, the other gaussian
          const double gamma = (*iter)->finite_nucleus() ? (*iter)->atom_exponent() : (*titer)->atom_exponent();
          out += charge * erf(sqrt(gamma) * dist) / dist;
        }
#endif
      }
    }
  }
  return out;
}


void Molecule::print_atoms() const {
  cout << "  *** Geometry ***" << endl << endl;
  for (auto i : atoms_) i->print();
  cout << endl;
}


array<double,3> Molecule::charge_center() const {
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


array<double,6> Molecule::quadrupole() const {
  array<double,6> out;
  array<double,3> c = charge_center();
  for (auto& i : atoms_) {
    out[0] += i->atom_charge() * pow(i->position(0) - c[0], 2);
    out[1] += i->atom_charge() * (i->position(0) - c[0]) * (i->position(1) - c[1]);
    out[2] += i->atom_charge() * pow(i->position(1) - c[1], 2);
    out[3] += i->atom_charge() * (i->position(0) - c[0]) * (i->position(2) - c[2]);
    out[4] += i->atom_charge() * (i->position(1) - c[1]) * (i->position(2) - c[2]);
    out[5] += i->atom_charge() * pow(i->position(2) - c[2], 2);
  }
  return out;
}


bool Molecule::has_finite_nucleus() const {
  return any_of(atoms_.begin(), atoms_.end(), [](shared_ptr<const Atom> a) { return a->finite_nucleus(); });
}


shared_ptr<const XYZFile> Molecule::xyz() const {
  auto out = make_shared<XYZFile>(natom());
  int iat = 0;
  for (auto& i : atoms_) {
    out->element(0, iat) = i->position(0);
    out->element(1, iat) = i->position(1);
    out->element(2, iat) = i->position(2);
    ++iat;
  }
  return out;
}


void Molecule::common_init1() {
  lmax_ = 0;
  aux_lmax_ = 0;
  nbasis_ = 0;
  naux_ = 0;
  nele_ = 0;
  nfrc_ = 0;

  for (auto& catom : atoms_) {
    nele_ += atommap.atom_number(catom->name());
    if (catom->use_ecp_basis()) nele_ -= catom->ecp_parameters()->ecp_ncore(); // valence electrons only

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


void Molecule::merge_obs_aux() {
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


int Molecule::num_count_ncore_only() const {
  int out = 0;
  for (auto& it : atoms_) {
    if (it->atom_number() >= 2) out += 2;
    if (it->atom_number() >= 10) out += 8;
    if (it->atom_number() >= 18) out += 8;
    if (it->atom_number() >= 36) out += 18;
    if (it->atom_number() >= 54) out += 18;
    if (it->atom_number() >= 86) out += 32;
    if (it->atom_number() >= 118) out += 32;
    if (it->atom_number() >  118) throw logic_error("Molecule::num_count_ncore_only() thinks you are using an atom with Z > 118...");
  }
  return out;
}


bool Molecule::operator==(const Molecule& o) const {
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

array<shared_ptr<const Matrix>,3> Molecule::compute_internal_coordinate(shared_ptr<const Matrix> prev,
    vector<shared_ptr<const OptExpBonds>> explicit_bond, bool negative_hessian,  bool verbose) const {
  if (verbose)
    cout << "    o Connectivitiy analysis" << endl;

  // list of primitive internals
  array<vector<vector<int>>,4> prims;
  vector<vector<double>> out;
  const size_t size = natom()*3;

  list<shared_ptr<Node>> nodes;
  int n = 0;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i, ++n) {
    if ((*i)->dummy()) throw runtime_error("haven't thought about internal coordinate with dummy atoms (or gradient in general)");
    nodes.push_back(make_shared<Node>(*i, n));
  }

  vector<double> hessprim;
  hessprim.reserve(natom()*3 * 10);

  // first pick up bonds
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    const double radiusi = (*i)->atom()->cov_radius();
    auto j = i;
    for (++j ; j != nodes.end(); ++j) {
      const double radiusj = (*j)->atom()->cov_radius();

      bool expbond = false;
      if (explicit_bond.size())
        for (auto e : explicit_bond)
          if ((e->pair(0) == (*i)->num()) && (e->pair(1) == (*j)->num())) expbond = true;

      if (((*i)->atom()->distance((*j)->atom()) < (radiusi+radiusj)*1.3) || expbond) {
        (*i)->add_connected(*j);
        (*j)->add_connected(*i);
        if (verbose)
        cout << "       bond:  " << setw(6) << (*i)->num() << setw(6) << (*j)->num() << "     bond length" <<
                                    setw(10) << setprecision(4) << (*i)->atom()->distance((*j)->atom()) << " bohr" << endl;

        // see IJQC 106, 2536 (2006)
        const double modelhess = 0.35 * adf_rho(*i, *j);
        hessprim.push_back(modelhess);

        Quatern<double> ip = (*i)->atom()->position();
        Quatern<double> jp = (*j)->atom()->position();

        jp -= ip;  // jp is a vector from i to j
        jp.normalize();
        vector<double> current(size);
        const double fac = adf_rho(*i, *j);

        // use baker, temporarily
        current[3*(*i)->num()+0] =  jp[1]*fac;
        current[3*(*i)->num()+1] =  jp[2]*fac;
        current[3*(*i)->num()+2] =  jp[3]*fac;
        current[3*(*j)->num()+0] = -jp[1]*fac;
        current[3*(*j)->num()+1] = -jp[2]*fac;
        current[3*(*j)->num()+2] = -jp[3]*fac;
        out.push_back(current);

        vector<int> prim(2);
        prim[0] = (*i)->num();
        prim[1] = (*j)->num();
        prims[0].push_back(prim);
      }
    }
  }

  if (!molecule_details::is_one_molecule(nodes))
    throw runtime_error("Internal coordinate transformation currently requires that we have only one molecule -- consider using \"explicitbond\".");

  // then bond angles A-O-B (A<B)
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    auto j = i;
    for (++j; j != nodes.end(); ++j) {
      std::set<shared_ptr<Node>> center = (*i)->common_center(*j);
      for (auto c = center.begin(); c != center.end(); ++c) {
        const double theta = (*c)->atom()->angle((*i)->atom(), (*j)->atom());
#if 0
        cout << "       angle: " << setw(6) << (*c)->num() << setw(6) << (*i)->num() << setw(6) << (*j)->num() <<
                "     angle" << setw(10) << setprecision(4) << theta << " deg" << endl;
#endif
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
        Quatern<double> st1 = (e21 * cos(rad) - e23) / (r21 * sin(rad));
        Quatern<double> st3 = (e23 * cos(rad) - e21) / (r23 * sin(rad));
        Quatern<double> st2 = (st1 + st3) * (-1.0);
        vector<double> current(size);
        // see IJQC 106, 2536 (2006)
        const double modelhess = 0.15 * adf_rho(*i, *c) * adf_rho(*c, *j);
        hessprim.push_back(modelhess);
        const double fval = 0.12;
        const double fac = sqrt(adf_rho(*i, *c) * adf_rho(*c, *j)) * (fval + (1-fval)*sin(rad));
        for (int ic = 0; ic != 3; ++ic) {
          current[3*(*i)->num() + ic] = st1[ic+1]*fac;
          current[3*(*j)->num() + ic] = st3[ic+1]*fac;
          current[3*(*c)->num() + ic] = st2[ic+1]*fac;
        }
        out.push_back(current);

        vector<int> prim(3);
        prim[0] = (*i)->num();
        prim[1] = (*c)->num();
        prim[2] = (*j)->num();
        prims[1].push_back(prim);
      }
    }
  }

  // proper dihedral angle (i-c-j-k) or improper dihedral angle (c-i,j,k)
  for (auto i = nodes.begin(); i != nodes.end(); ++i) {
    for (auto j = nodes.begin(); j != nodes.end(); ++j) {
      if (*i == *j) continue;
      std::set<shared_ptr<Node>> center = (*i)->common_center(*j);
      for (auto k = nodes.begin(); k != nodes.end(); ++k) {
        for (auto c = center.begin(); c != center.end(); ++c) {
          bool improper;
          if (((*k)->connected_with(*j) && (*j)->connected_with(*c) && (*c)->connected_with(*i))) {
            if (*c > *j) continue;
            if (*c == *k || *k == *i) continue;
            if (*j == *k) continue;
            improper = false;
          } else if (((*k)->connected_with(*c) && (*j)->connected_with(*c) && (*i)->connected_with(*c))) {
            if (*j == *k || *k == *i) continue;
            if (!(*i < *j && *j < *k)) continue;
            improper = true;
          } else continue;
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
          const double tabc = atan2(rotabc.norm(), -eab.dot_product(ebc));
          const double tbcd = atan2(rotbcd.norm(), -edc.dot_product(ecb));

          Quatern<double> sa = (eab * ebc) / (-rab*::pow(::sin(tabc), 2.0));
          Quatern<double> sd = (edc * ecb) / (-rcd*::pow(::sin(tbcd), 2.0));
          Quatern<double> sb = (eab * ebc) * ((rbc-rab*::cos(tabc)) / (rab*rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * (::cos(tbcd) / (rbc*::pow(::sin(tbcd), 2.0)));
          Quatern<double> sc = (eab * ebc) * (::cos(tabc) / (rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * ((rbc-rcd*::cos(tbcd)) / (rcd*rbc*::pow(::sin(tbcd), 2.0)));
          vector<double> current(size);
          // see IJQC 106, 2536 (2006)
          double modelhess;
          if (improper) modelhess = 0.005 * adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*c, *k);
          else          modelhess = 0.005 * adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*j, *k);
          hessprim.push_back(modelhess);
          const double  theta0 = (*c)->atom()->angle((*i)->atom(), (*j)->atom()) / rad2deg__;
          double theta1;
          if (improper) theta1 = (*c)->atom()->angle((*i)->atom(), (*k)->atom()) / rad2deg__;
          else          theta1 = (*j)->atom()->angle((*c)->atom(), (*k)->atom()) / rad2deg__;
          const double fval = 0.12;
          double fac;
          if (improper) fac = pow(adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*c, *k), 1.0/3.0) * (fval + (1-fval)*sin(theta0)) * (fval + (1-fval)*sin(theta1));
          else          fac = pow(adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*j, *k), 1.0/3.0) * (fval + (1-fval)*sin(theta0)) * (fval + (1-fval)*sin(theta1));
          for (int ic = 0; ic != 3; ++ic) {
            current[3*(*i)->num() + ic] = sa[ic+1]*fac;
            current[3*(*c)->num() + ic] = sb[ic+1]*fac;
            current[3*(*j)->num() + ic] = sc[ic+1]*fac;
            current[3*(*k)->num() + ic] = sd[ic+1]*fac;
            assert(fabs(sa[ic+1]+sb[ic+1]+sc[ic+1]+sd[ic+1]) < 1.0e-8);
          }
          out.push_back(current);
          vector<int> prim(4);
          prim[0] = (*i)->num();
          prim[1] = (*c)->num();
          prim[2] = (*j)->num();
          prim[3] = (*k)->num();
          if (improper)  prims[3].push_back(prim);
          else           prims[2].push_back(prim);
        }
      }
    }
  }

  // debug output
  const int primsize = out.size();
  const int cartsize = 3*natom();

  vector<vector<double>> cvec;
  // bdag = B^T
  Matrix bdag(cartsize, primsize);
  double* biter = bdag.data();
  for (auto i = out.begin(); i != out.end(); ++i, biter += cartsize)
    copy(i->begin(), i->end(), biter);

  // TODO this is needed but I don't know why..
  bdag.broadcast();

  Matrix bb = bdag % bdag;
  VectorB eig(primsize);
  bb.diagonalize(eig);

  // make them consistent if parallel and not using ScaLapack
#ifndef HAVE_SCALAPACK
  bb.broadcast();
  mpi__->broadcast(eig.data(), primsize, 0);
#endif


  int ninternal = max(cartsize-6,1);
  for (int i = primsize-ninternal; i != primsize; ++i) {
    if (eig(i) < 1.0e-10)
      cout << "       ** caution **  small eigenvalue " << eig(i) << endl;
  }
  if (verbose)
    cout << "      Nonredundant internal coordinate generated (dim = " << ninternal << ")" << endl;

  auto umat = make_shared<Matrix>(bb.slice(primsize-ninternal,primsize));
  // B = U^T * B^prim, B^T = (B^prim)^T * U
  auto bnew = make_shared<Matrix>(bdag * *umat);

  // form (B^+)^-1 = (BB^+)^-1 B = Lambda^-1 B
  auto bdmnew = make_shared<Matrix>(cartsize, ninternal);
  for (int i = 0; i != ninternal; ++i)
    for (int j = 0; j != cartsize; ++j)
      bdmnew->element(j,i) = bnew->element(j,i) / eig[i+primsize-ninternal];

  // compute hessian
  auto hessp = make_shared<Matrix>(primsize, primsize);
  for (int i = 0; i != primsize; ++i)
    hessp->element(i,i) = hessprim[i];

  if (negative_hessian) {
    for (int i = 0; i != primsize; ++i)
      for (int j = i+1; j != primsize; ++j) {
        hessp->element(i,j) = -sqrt(2.0 * hessprim[i] * hessprim[j]);
        hessp->element(j,i) = -sqrt(2.0 * hessprim[i] * hessprim[j]);
      }
  }

  auto hess = make_shared<Matrix>(*umat % *hessp * *umat);
  if (!negative_hessian) {
    hess->sqrt();
    *bnew = *bnew * *hess;
    hess->inverse();
    *bdmnew = *bdmnew * *hess;
  }
  auto hessout = make_shared<Matrix>(*hess);

  if (prev) {
    // internal--internal matrix
    Matrix approx1 = *prev % *bdmnew;
    assert(approx1.ndim() == ninternal && approx1.mdim() == ninternal);
    *bnew = *prev;
    try {
      approx1.inverse();
    } catch (const runtime_error& error) {
      throw runtime_error("It seems that the geometry has changed substantially. Start over the optimization.");
    }
    *bdmnew *= approx1;
  }

  // make them consistent
  bnew->broadcast();
  bdmnew->broadcast();

  return array<shared_ptr<const Matrix>,3>{{bnew, bdmnew, hessout}};
}

tuple<vector<vector<int>>,array<shared_ptr<const Matrix>,5>> Molecule::compute_redundant_coordinate(vector<vector<int>> prev_bond) const {
  vector<vector<double>> out;

  vector<vector<int>> bondlist;
  vector<double> val;
  const size_t size = natom()*3;

  list<shared_ptr<Node>> nodes;
  int n = 0;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i, ++n) {
    if ((*i)->dummy()) throw runtime_error("haven't thought about internal coordinate with dummy atoms (or gradient in general)");
    nodes.push_back(make_shared<Node>(*i, n));
  }

  vector<double> hessprim;
  hessprim.reserve(natom()*3 * 10);

  if (prev_bond.size() == 0) {
    // first pick up bonds. this is done only once
    cout << "    o Connectivitiy analysis" << endl;
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
      const double radiusi = (*i)->atom()->cov_radius();
      auto j = i;
      for (++j ; j != nodes.end(); ++j) {
        const double radiusj = (*j)->atom()->cov_radius();

        if ((*i)->atom()->distance((*j)->atom()) < (radiusi+radiusj)*1.3) {
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
          hessprim.push_back(0.50);

          const double length = (*i)->atom()->distance((*j)->atom()) * -1.0;
          val.push_back(length);

          vector<int> bondij(3);
          bondij[0] = 0;
          bondij[1] = (*i)->num();
          bondij[2] = (*j)->num();
          bondlist.push_back(bondij);
        }
      }
    }

    if (!molecule_details::is_one_molecule(nodes))
      throw runtime_error("Internal coordinate transformation currently requires that we have only one molecule -- consider using \"explicitbond\".");

    // then bond angles A-O-B (A<B)
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
      auto j = i;
      for (++j; j != nodes.end(); ++j) {
        std::set<shared_ptr<Node>> center = (*i)->common_center(*j);
        for (auto c = center.begin(); c != center.end(); ++c) {
          const double theta = (*c)->atom()->angle((*i)->atom(), (*j)->atom());
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
          Quatern<double> st1 = (e21 * cos(rad) - e23) / (r21 * sin(rad));
          Quatern<double> st3 = (e23 * cos(rad) - e21) / (r23 * sin(rad));
          Quatern<double> st2 = (st1 + st3) * (-1.0);
          vector<double> current(size);

          for (int ic = 0; ic != 3; ++ic) {
            current[3*(*i)->num() + ic] = st1[ic+1];
            current[3*(*j)->num() + ic] = st3[ic+1];
            current[3*(*c)->num() + ic] = st2[ic+1];
          }
          out.push_back(current);
          hessprim.push_back(0.20);

          val.push_back(rad);

          vector<int> bondij(4);
          bondij[0] = 1;
          bondij[1] = (*i)->num();
          bondij[2] = (*j)->num();
          bondij[3] = (*c)->num();
          bondlist.push_back(bondij);
        }
      }
    }

    // then dihedral angle (i-c-j-k)
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
      for (auto j = nodes.begin(); j != nodes.end(); ++j) {
        if (*i == *j) continue;
        std::set<shared_ptr<Node>> center = (*i)->common_center(*j);
        for (auto k = nodes.begin(); k != nodes.end(); ++k) {
          for (auto c = center.begin(); c != center.end(); ++c) {
            bool improper;
            if (((*k)->connected_with(*j) && (*j)->connected_with(*c) && (*c)->connected_with(*i))) {
              if (*c > *j) continue;
              if (*c == *k || *k == *i) continue;
              if (*j == *k) continue;
              improper = false;
            } else if (((*k)->connected_with(*c) && (*j)->connected_with(*c) && (*i)->connected_with(*c))) {
              if (*j == *k || *k == *i) continue;
              if (!(*i < *j && *j < *k)) continue;
              improper = true;
            } else continue;
            double dihed = (*c)->atom()->dihedral_angle((*i)->atom(), (*j)->atom(), (*k)->atom())/rad2deg__;
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
            const double tabc = atan2(rotabc.norm(), -eab.dot_product(ebc));
            const double tbcd = atan2(rotbcd.norm(), -edc.dot_product(ecb));

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
              assert(fabs(sa[ic+1]+sb[ic+1]+sc[ic+1]+sd[ic+1]) < 1.0e-8);
            }
            out.push_back(current);
            hessprim.push_back(0.1);

            val.push_back(dihed);

            vector<int> bondij(6);
            bondij[0] = 2;
            bondij[1] = (*i)->num();
            bondij[2] = (*c)->num();
            bondij[3] = (*j)->num();
            bondij[4] = (*k)->num();
            bondij[5] = improper ? 1 : 0;
            bondlist.push_back(bondij);
          }
        }
      }
    }
    cout << "      Redundant internal coordinate generated (dim = " << bondlist.size() << ")" << endl;

  } else {

    // Now I have original connectivity matrix. Never update the bondlist when redundant

    for (int i = 0; i != prev_bond.size(); ++i) {
      if (prev_bond[i][0] == 0) {
        const int iatom = prev_bond[i][1];
        const int jatom = prev_bond[i][2];
        Quatern<double> ip = atoms_[iatom]->position();
        Quatern<double> jp = atoms_[jatom]->position();
        jp -= ip;  // jp is a vector from i to j
        jp.normalize();
        vector<double> current(size);
        current[3*iatom+0] =  jp[1];
        current[3*iatom+1] =  jp[2];
        current[3*iatom+2] =  jp[3];
        current[3*jatom+0] = -jp[1];
        current[3*jatom+1] = -jp[2];
        current[3*jatom+2] = -jp[3];
        out.push_back(current);
        hessprim.push_back(0.50);

        const double length = atoms_[iatom]->distance(atoms_[jatom]) * -1.0;
        val.push_back(length);
      } else if (prev_bond[i][0] == 1) {
        const int iatom = prev_bond[i][1];
        const int jatom = prev_bond[i][2];
        const int catom = prev_bond[i][3];

        const double theta = atoms_[catom]->angle(atoms_[iatom], atoms_[jatom]);
        Quatern<double> op = atoms_[catom]->position();
        Quatern<double> ap = atoms_[iatom]->position();
        Quatern<double> bp = atoms_[jatom]->position();
        Quatern<double> e21 = ap - op;
        Quatern<double> e23 = bp - op;
        const double r21 = e21.norm();
        const double r23 = e23.norm();
        e21.normalize();
        e23.normalize();
        const double rad = theta/rad2deg__;
        Quatern<double> st1 = (e21 * cos(rad) - e23) / (r21 * sin(rad));
        Quatern<double> st3 = (e23 * cos(rad) - e21) / (r23 * sin(rad));
        Quatern<double> st2 = (st1 + st3) * (-1.0);
        vector<double> current(size);

        for (int ic = 0; ic != 3; ++ic) {
          current[3*iatom + ic] = st1[ic+1];
          current[3*jatom + ic] = st3[ic+1];
          current[3*catom + ic] = st2[ic+1];
        }
        out.push_back(current);
        hessprim.push_back(0.20);

        val.push_back(rad);
      } else if (prev_bond[i][0] == 2) {
        const int iatom = prev_bond[i][1];
        const int catom = prev_bond[i][2];
        const int jatom = prev_bond[i][3];
        const int katom = prev_bond[i][4];

        double dihed = atoms_[catom]->dihedral_angle(atoms_[iatom], atoms_[jatom], atoms_[katom])/rad2deg__;
        Quatern<double> ap = atoms_[iatom]->position();
        Quatern<double> bp = atoms_[catom]->position();
        Quatern<double> cp = atoms_[jatom]->position();
        Quatern<double> dp = atoms_[katom]->position();
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
        const double tabc = atan2(rotabc.norm(), -eab.dot_product(ebc));
        const double tbcd = atan2(rotbcd.norm(), -edc.dot_product(ecb));

        Quatern<double> sa = (eab * ebc) / (-rab*::pow(::sin(tabc), 2.0));
        Quatern<double> sd = (edc * ecb) / (-rcd*::pow(::sin(tbcd), 2.0));
        Quatern<double> sb = (eab * ebc) * ((rbc-rab*::cos(tabc)) / (rab*rbc*::pow(::sin(tabc), 2.0)))
                           + (edc * ecb) * (::cos(tbcd) / (rbc*::pow(::sin(tbcd), 2.0)));
        Quatern<double> sc = (eab * ebc) * (::cos(tabc) / (rbc*::pow(::sin(tabc), 2.0)))
                           + (edc * ecb) * ((rbc-rcd*::cos(tbcd)) / (rcd*rbc*::pow(::sin(tbcd), 2.0)));
        vector<double> current(size);

        for (int ic = 0; ic != 3; ++ic) {
          current[3*iatom + ic] = sa[ic+1];
          current[3*catom + ic] = sb[ic+1];
          current[3*jatom + ic] = sc[ic+1];
          current[3*katom + ic] = sd[ic+1];
        }
        out.push_back(current);
        hessprim.push_back(0.10);

        val.push_back(dihed);
      }
    }
    bondlist = prev_bond;
  }

  // debug output
  const int primsize = out.size();
  const int cartsize = 3*natom();

  auto inithess = make_shared<Matrix>(primsize, primsize);
  for (int i = 0; i != primsize; ++i) {
    inithess->element(i, i) = hessprim[i];
  }

  Matrix minv(cartsize, cartsize);
  minv.zero();
  for (int i = 0; i != natom(); ++i) {
    minv(i*3+0,i*3+0) = 1.0;
    minv(i*3+1,i*3+1) = 1.0;
    minv(i*3+2,i*3+2) = 1.0;
  }

  // By convention, bdag is B^+ here
  Matrix bdag(cartsize, primsize);
  double* biter = bdag.data();
  for (auto i = out.begin(); i != out.end(); ++i, biter += cartsize)
    copy(i->begin(), i->end(), biter);

  bdag.broadcast();

  auto bmat = make_shared<Matrix>(*bdag.transpose());
  bmat->broadcast();


  Matrix g = bdag % minv * bdag;
  shared_ptr<Matrix> g2 = g.copy();
  VectorB geig(primsize);
  // (-1.0) is a trick to make it (K L), not (L K) (the results are essentially the same, but just following the literature...)
  g.scale(-1.0);
  g.diagonalize(geig);

#ifndef HAVE_SCALAPACK
  g.broadcast();
  mpi__->broadcast(geig.data(), primsize, 0);
#endif

  // pseudoinvert g
  Matrix ml(primsize, primsize);
  ml.zero();
  auto values = make_shared<Matrix>(3, primsize/3+1);

  for (int i = 0; i != primsize; ++i) {
    values->element(i%3,i/3) = val[i];
    geig(i) *= -1.0;
    if (fabs(geig(i)) < 1.0e-7) { geig(i) = 0.0; ml(i,i) = 0.0;
      for (int j = 0; j != primsize; ++j)
        g(j,i) = 0.0;
    }
    else { ml(i,i) = 1.0 / geig(i); }
  }
  Matrix ginv = g * ml ^ g;
  auto ubgnew = make_shared<Matrix>(minv * bdag * ginv);
  auto gg = make_shared<Matrix>(*g2 * ginv);

  ubgnew->broadcast();

  return make_tuple(bondlist,array<shared_ptr<const Matrix>,5>{{bmat, ubgnew, values, inithess, gg}});
}


const vector<shared_ptr<const Molecule>> Molecule::split_atoms(const int max_atoms) const {
  vector<shared_ptr<const Molecule>> out = {};
  const int res = natom() % max_atoms;
  const int nfullset = natom() / max_atoms;
  const int nset = nfullset + (res != 0 ? 1 : 0);
  assert((nset == nfullset && res == 0) || (nset == nfullset+1 && res > 0));
  assert(nfullset*max_atoms + res == natom());
  out.resize(nset);

  // Do not bother copying aux_atom data, since they are not needed for NAI
  const vector<shared_ptr<const Atom>> empty_aux_atoms = {};

  for (int i=0; i!=nfullset; ++i) {
    vector<shared_ptr<const Atom>> current_atoms(&atoms_[i*max_atoms], &atoms_[(i+1)*max_atoms]);
    auto current_mol = make_shared<const Molecule>(current_atoms, empty_aux_atoms);
    out[i] = current_mol;
  }
  if (res != 0) {
    vector<shared_ptr<const Atom>> current_atoms(&atoms_[nfullset*max_atoms], &atoms_[nfullset*max_atoms+res]);
    auto current_mol = make_shared<const Molecule>(current_atoms, empty_aux_atoms);
    out[nfullset] = current_mol;
  }
  return out;
}


shared_ptr<Molecule> Molecule::uncontract() const {
  vector<shared_ptr<const Atom>> atom;
  for (auto& i : atoms_)
    atom.push_back(i->uncontract()->relativistic());

  auto mol = make_shared<Molecule>(atom, aux_atoms_);
  mol->spherical_ = spherical_;
  mol->aux_merged_ = aux_merged_;
  mol->basisfile_ = basisfile_;
  mol->auxfile_ = auxfile_;
  mol->nuclear_repulsion_ = nuclear_repulsion_;
  mol->external_ = external_; 
  mol->magnetic_field_ = magnetic_field_; 
  mol->skip_self_interaction_ = skip_self_interaction_; 
  mol->cap_ = cap_;

  mol->common_init1();

  return mol;
}
