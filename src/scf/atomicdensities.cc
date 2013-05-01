//
// BAGEL - Parallel electron correlation program.
// Filename: atomicdensities.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/scf/atomicdensities.h>
#include <src/util/atommap.h>
#include <src/util/diis.h>
#include <src/scf/rohf.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap_;

AtomicDensities::AtomicDensities(std::shared_ptr<const Geometry> g) : Matrix(g->nbasis(), g->nbasis()), geom_(g) {
  // first make a list of unique atoms
  const string basis = geom_->basisfile();
  const string dfbasis = geom_->auxfile();
  map<string, shared_ptr<const Matrix>> atoms;
  unique_ptr<double[]> eig(new double[geom_->nbasis()]);

  int offset = 0;
  for (auto& i : geom_->atoms()) {
    if (i->dummy()) continue;
    if (atoms.find(i->name()) == atoms.end()) {
      // dummy buffer to suppress the output
      stringstream ss;
      std::streambuf* cout_orig = cout.rdbuf();
      cout.rdbuf(ss.rdbuf());

      shared_ptr<const Atom> atom(new Atom(i->spherical(), i->name(), {{0.0,0.0,0.0}}, basis));

      boost::property_tree::ptree geomop;
      geomop.put("basis", basis);
      geomop.put("df_basis", dfbasis.empty() ? basis : dfbasis);
      shared_ptr<const Geometry> ga(new Geometry({atom}, geomop));
      atoms.insert(make_pair(i->name(), compute_atomic(ga)));

      // restore cout
      cout.rdbuf(cout_orig);
    }
    auto iter = atoms.find(i->name());
    assert(iter != atoms.end());
    copy_block(offset, offset, i->nbasis(), i->nbasis(), iter->second->data());
    offset += i->nbasis();
  }

}


shared_ptr<const Matrix> AtomicDensities::compute_atomic(shared_ptr<const Geometry> ga) const {
  // this thing does not work with cartesian basis
  assert(ga->spherical());
  // first the number of s, p, d, f orbitals
  shared_ptr<const Atom> atom = ga->atoms().front();
  vector<int> num(5);
  for (auto& i : atom->shells()) {
    const int ang = i->angular_number();
    if (ang < 4) num[ang] += i->nbasis();
  }
  num[4] = ga->nbasis() - accumulate(&num[0], &num[4], 0);

  // overlap matrix
  shared_ptr<Overlap> overlap(new Overlap(ga));
  // block diagonal structure is (should be)  maintained
  TildeX tildex(overlap, 1.0e-5);
  shared_ptr<Hcore> hcore(new Hcore(ga));
  DIIS<Matrix> diis(5);
  shared_ptr<const Matrix> coeff;
  unique_ptr<double[]> eig(new double[ga->nbasis()]);
  {
    Matrix ints = tildex % *hcore * tildex;
    ints = *ints.diagonalize_blocks(eig.get(), num);
    coeff = shared_ptr<const Matrix>(new Matrix(tildex * ints));
  }

  tuple<int,int,int,int> nclosed = atommap_.num_closed(ga->atoms().front()->name());
  tuple<int,int,int,int> nopen   = atommap_.num_open(ga->atoms().front()->name());
  int nclo[4] = {get<0>(nclosed)/2, get<1>(nclosed)/2, get<2>(nclosed)/2, get<3>(nclosed)/2};
  int nope[4] = {get<0>(nopen),     get<1>(nopen),     get<2>(nopen),     get<3>(nopen)};
  const int sclosed = get<0>(nclosed)+get<1>(nclosed)+get<2>(nclosed)+get<3>(nclosed);
  const int sopen   = get<0>(nopen)  +get<1>(nopen)  +get<2>(nopen)  +get<3>(nopen);
  if (sclosed+sopen != ga->nele())
    throw logic_error("Inconsistent nclosed and nopen. See AtomMap");

  shared_ptr<Matrix> ocoeff(new Matrix(ga->nbasis(), ga->nbasis()));
  shared_ptr<Matrix> vcoeff(new Matrix(ga->nbasis(), ga->nbasis()));

  vector<int> orb{1,3,5,7};

  for (int i = 0, j = 0, k = 0, l = 0; i != 4; ++i) {
    for (int jj = 0; jj != nclo[i]; ++jj)
      copy_n(coeff->element_ptr(0,k+jj), ga->nbasis(), ocoeff->element_ptr(0,j+jj));
    if (nope[i]) {
      const int n = orb[i];
      const double occupation = std::sqrt(static_cast<double>(nope[i])/n);
      if (nclo[i]+n > num[i]) throw runtime_error("The basis set is smaller than the minimal basis set for " + atom->name() + ". Please check.");
      for (int jj = nclo[i]; jj != nclo[i]+n; ++jj) {
        daxpy_(ga->nbasis(), occupation, coeff->element_ptr(0,k+jj), 1, vcoeff->element_ptr(0,l++),1);
      }
    }
    j += nclo[i];
    k += num[i];
  }

  shared_ptr<Matrix> vden(new Matrix(*vcoeff ^ *vcoeff));

  int iter = 0;
  const int maxiter = 100;
  double prev_energy = 0.0;
  for (; iter != maxiter; ++iter) {
    shared_ptr<Matrix> fock = sclosed ? shared_ptr<Matrix>(new Fock<1>(ga, hcore, std::shared_ptr<const Matrix>(), ocoeff->slice(0,sclosed/2), true)) : hcore->copy();
    shared_ptr<Matrix> fock2(new Fock<1>(ga, hcore, vden, std::shared_ptr<const Matrix>(), false, 0.0));
    *fock += *fock2 - *hcore;

    shared_ptr<const Matrix> aodensity(new Matrix((*ocoeff^*ocoeff)*2.0 + *vden));
    const double energy = ((*hcore+*fock) * *aodensity).trace()*0.5;
    cout << setprecision(10) << energy << endl;

    shared_ptr<const Matrix> residual(new Matrix(*fock**aodensity**overlap - *overlap**aodensity**fock));
    if (residual->rms() < 1.0e-2 || fabs(energy-prev_energy) < 1.0e-4) break;
    prev_energy = energy;
    fock = diis.extrapolate(make_pair(fock, residual));

    Matrix ints = tildex % *fock * tildex;
    ints = *ints.diagonalize_blocks(eig.get(), num);

    coeff  = shared_ptr<const Matrix>(new Matrix(tildex * ints));
    ocoeff = shared_ptr<Matrix>(new Matrix(ga->nbasis(), sclosed));
    vcoeff = shared_ptr<Matrix>(new Matrix(ga->nbasis(), ga->nbasis()));
    for (int i = 0, j = 0, k = 0, l = 0; i != 4; ++i) {
      for (int jj = 0; jj != nclo[i]; ++jj)
        copy_n(coeff->element_ptr(0,k+jj), ga->nbasis(), ocoeff->element_ptr(0,j+jj));
      if (nope[i]) {
        const int n = orb[i];
        const double occupation = std::sqrt(static_cast<double>(nope[i])/n);
        for (int jj = nclo[i]; jj != nclo[i]+n; ++jj) {
          daxpy_(ga->nbasis(), occupation, coeff->element_ptr(0,k+jj), 1, vcoeff->element_ptr(0,l++),1);
        }
      }
      j += nclo[i];
      k += num[i];
    }
    vden = shared_ptr<Matrix>(new Matrix(*vcoeff ^ *vcoeff));
  }
  if (iter == maxiter) throw runtime_error("spin-averaged atomic HF did not converge");

  shared_ptr<Matrix> out = vden;
  if (sclosed) {
    shared_ptr<const Coeff> c(new Coeff(*ocoeff));
    *out += *c->form_density_rhf(ocoeff->mdim());
  }
  return out;
}
