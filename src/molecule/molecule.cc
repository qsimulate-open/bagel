//
// BAGEL - Parallel electron correlation program.
// Filename: molecule.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/wfn/geometry_connect.h>
#include <src/molecule/molecule.h>
#include <src/util/atommap.h>
#include <src/math/quatern.h>

using namespace std;
using namespace bagel;
using geometry_details::Node;
using geometry_details::adf_rho;

const static AtomMap atommap_;


BOOST_CLASS_EXPORT_IMPLEMENT(Molecule)

double Molecule::compute_nuclear_repulsion() {
  double out = 0.0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    double c = (*iter)->atom_charge();
    if ((*iter)->use_ecp_basis()) c -= (*iter)->ecp_parameters()->ecp_ncore();
    for (auto titer = iter + 1; titer != atoms_.end(); ++titer) {
      // nuclear repulsion between dummy atoms are not computed here (as in Molpro)
      if (!(*iter)->dummy() || !(*titer)->dummy()) {
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
          const double gamma0 = (*iter)->finite_nucleus();
          const double gamma1 = (*titer)->finite_nucleus();
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
  cout << "  Symmetry: " << symmetry() << endl;
  cout << endl;
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
    out[2] += i->atom_charge() * (i->position(0) - c[0]) * (i->position(2) - c[2]);
    out[3] += i->atom_charge() * pow(i->position(1) - c[1], 2);
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
    nele_ += atommap_.atom_number(catom->name());
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


int Molecule::num_count_ncore_only() const {
  int out = 0;
  for (auto& it : atoms_) {
    if (it->atom_number() >= 2) out += 2;
    if (it->atom_number() >= 10) out += 8;
    if (it->atom_number() >= 18) out += 8;
    if (it->atom_number() > 36) throw logic_error("needs to modify Molecule::count_num_ncore for atoms beyond Kr"); // TODO
  }
  return out;
}


int Molecule::num_count_full_valence_nocc() const {
  int out = 0;
  for (auto& it : atoms_) {
    if (it->atom_number() < 2) out += 1;
    if (it->atom_number() >= 2 && it->atom_number() <= 10) out += 5;
    if (it->atom_number() > 10) throw logic_error("needs to modify Molecule::num_count_full_valence_nocc for atoms beyond Ne"); // TODO
  }
  return out;
};


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


array<shared_ptr<const Matrix>,2> Molecule::compute_internal_coordinate(shared_ptr<const Matrix> prev) const {
  cout << "    o Connectivitiy analysis" << endl;

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

      if ((*i)->atom()->distance((*j)->atom()) < (radiusi+radiusj)*1.3) {
        (*i)->add_connected(*j);
        (*j)->add_connected(*i);
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
        current[3*(*i)->num()+0] =  jp[1]*fac;
        current[3*(*i)->num()+1] =  jp[2]*fac;
        current[3*(*i)->num()+2] =  jp[3]*fac;
        current[3*(*j)->num()+0] = -jp[1]*fac;
        current[3*(*j)->num()+1] = -jp[2]*fac;
        current[3*(*j)->num()+2] = -jp[3]*fac;
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
        Quatern<double> st1 = (e21 * ::cos(rad) - e23) / (r21 * ::sin(rad));
        Quatern<double> st3 = (e23 * ::cos(rad) - e21) / (r23 * ::sin(rad));
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
#if 0
          cout << "    dihedral: " << setw(6) << (*i)->num() << setw(6) << (*c)->num() << setw(6) << (*j)->num() << setw(6) << (*k)->num() <<
                  "     angle" << setw(10) << setprecision(4) << (*c)->atom()->dihedral_angle((*i)->atom(), (*j)->atom(), (*k)->atom()) << " deg" << endl;
#endif
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
          const double tabc = ::atan2(rotabc.norm(), -eab.dot_product(ebc));
          const double tbcd = ::atan2(rotbcd.norm(), -edc.dot_product(ecb));

          Quatern<double> sa = (eab * ebc) / (-rab*::pow(::sin(tabc), 2.0));
          Quatern<double> sd = (edc * ecb) / (-rcd*::pow(::sin(tbcd), 2.0));
          Quatern<double> sb = (eab * ebc) * ((rbc-rab*::cos(tabc)) / (rab*rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * (::cos(tbcd) / (rbc*::pow(::sin(tbcd), 2.0)));
          Quatern<double> sc = (eab * ebc) * (::cos(tabc) / (rbc*::pow(::sin(tabc), 2.0)))
                             + (edc * ecb) * ((rbc-rcd*::cos(tbcd)) / (rcd*rbc*::pow(::sin(tbcd), 2.0)));
          vector<double> current(size);
          // see IJQC 106, 2536 (2006)
          const double modelhess = 0.005 * adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*j, *k);
          hessprim.push_back(modelhess);
          const double theta0 = (*c)->atom()->angle((*i)->atom(), (*j)->atom()) / rad2deg__;
          const double theta1 = (*j)->atom()->angle((*c)->atom(), (*k)->atom()) / rad2deg__;
          const double fval = 0.12;
          const double fac = pow(adf_rho(*i, *c) * adf_rho(*c, *j) * adf_rho(*j, *k), 1.0/3.0) * (fval + (1-fval)*sin(theta0)) * (fval + (1-fval)*sin(theta1));
          for (int ic = 0; ic != 3; ++ic) {
            current[3*(*i)->num() + ic] = sa[ic+1]*fac;
            current[3*(*c)->num() + ic] = sb[ic+1]*fac;
            current[3*(*j)->num() + ic] = sc[ic+1]*fac;
            current[3*(*k)->num() + ic] = sd[ic+1]*fac;
            assert(fabs(sa[ic+1]+sb[ic+1]+sc[ic+1]+sd[ic+1]) < 1.0e-8);
          }
          out.push_back(current);
        }
      }
    }
  }

  // debug output
  const int primsize = out.size();
  const int cartsize = 3*natom();

  Matrix bdag(cartsize, primsize);
  double* biter = bdag.data();
  for (auto i = out.begin(); i != out.end(); ++i, biter += cartsize)
    copy(i->begin(), i->end(), biter);

  // TODO this is needed but I don't know why..
  bdag.broadcast();

  Matrix bb = bdag % bdag * (-1.0);
  VectorB eig(primsize);
  bb.diagonalize(eig);

  // make them consistent if parallel and not using ScaLapack
#ifndef HAVE_SCALAPACK
  bb.broadcast();
  mpi__->broadcast(eig.data(), primsize, 0);
#endif

  int ninternal = max(cartsize-6,1);
  for (int i = 0; i != ninternal; ++i) {
    eig(i) *= -1.0;
    if (eig(i) < 1.0e-10)
      cout << "       ** caution **  small eigenvalue " << eig(i) << endl;
  }
  cout << "      Nonredundant internal coordinate generated (dim = " << ninternal << ")" << endl;

  // form B = U^+ Bprim
  Matrix bbslice = bb.slice(0,ninternal);
  auto bnew = make_shared<Matrix>(bdag * bbslice);

  // form (B^+)^-1 = (BB^+)^-1 B = Lambda^-1 B
  auto bdmnew = make_shared<Matrix>(cartsize, ninternal);
  for (int i = 0; i != ninternal; ++i)
    for (int j = 0; j != cartsize; ++j)
      bdmnew->element(j,i) = bnew->element(j,i) / eig[i];

  // compute hessian
  Matrix scale = bbslice;
  for (int i = 0; i != ninternal; ++i) {
    for (int j = 0; j != primsize; ++j) {
      scale(j,i) *= hessprim[j];
    }
  }
  Matrix hess = bbslice % scale;
  hess.sqrt();
  *bnew = *bnew * hess;
  hess.inverse();
  *bdmnew = *bdmnew * hess;

  // if this is not the first time, make sure that the change is minimum
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

  return array<shared_ptr<const Matrix>,2>{{bnew, bdmnew}};
}

