//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: atomicdensities.cc
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

#include <src/scf/atomicdensities.h>
#include <src/scf/hf/rohf.h>
#include <src/scf/hf/fock.h>
#include <src/util/atommap.h>
#include <src/util/math/diis.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap;

AtomicDensities::AtomicDensities(shared_ptr<const Geometry> g) : Matrix(g->nbasis(), g->nbasis()), geom_(g) {
  // first make a list of unique atoms
  const string defbasis = geom_->basisfile();
  map<pair<string,string>, shared_ptr<const Matrix>> atoms;

  int offset = 0;

  // read basis file
  shared_ptr<const PTree> bdata = PTree::read_basis(defbasis);

  vector<shared_ptr<const Atom>> aux_atoms;
  if (geom_->auxfile().empty()) {
     for (auto& a : geom_->atoms()) {
       auto aux_atom = make_shared<const Atom>(*a, a->spherical(), geom_->basisfile(), make_pair(geom_->basisfile(), bdata), nullptr);
       aux_atoms.push_back(aux_atom);
     }
  } else {
    aux_atoms = geom_->aux_atoms();
  }
  auto ai = aux_atoms.begin();
  for (auto& i : geom_->atoms()) {
    if (i->dummy()) { ++ai; continue; }
    if (atoms.find({i->name(),i->basis()}) == atoms.end()) {
      // dummy buffer to suppress the output
      stringstream ss;
      std::streambuf* cout_orig = cout.rdbuf();
      cout.rdbuf(ss.rdbuf());

      shared_ptr<PTree> geomop = make_shared<PTree>();
      const string basis = i->basis();
      geomop->put("basis", basis);
      const string dfbasis = (*ai)->basis();
      geomop->put("df_basis", !dfbasis.empty() ? dfbasis : basis);

      auto atom = make_shared<const Atom>(i->spherical(), i->name(), array<double,3>{{0.0,0.0,0.0}}, basis, make_pair(defbasis, bdata), nullptr);
      // TODO geometry makes aux atoms, which is ugly
      auto ga = make_shared<const Geometry>(vector<shared_ptr<const Atom>>{atom}, geomop);
      atoms.emplace(make_pair(i->name(),i->basis()), compute_atomic(ga));

      // restore cout
      cout.rdbuf(cout_orig);
    }
    auto iter = atoms.find({i->name(),i->basis()});
    assert(iter != atoms.end());
    copy_block(offset, offset, i->nbasis(), i->nbasis(), iter->second);
    offset += i->nbasis();

    ++ai;
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
  num[4] = ga->nbasis() - accumulate(num.data(), num.data()+4, 0);

  // overlap matrix
  auto overlap = make_shared<Overlap>(ga);
  // block diagonal structure is (should be)  maintained
  // instead of tildex, here we use inverse_half (because we will not be able to use diagonalize_blocks)
  auto tildex = overlap->copy();
  tildex->inverse_half(1.0e-5);
  auto hcore = make_shared<Hcore>(ga);
  DIIS<Matrix> diis(5);
  shared_ptr<const Matrix> coeff;
  VectorB eig(ga->nbasis());
  {
    Matrix ints = *tildex % *hcore * *tildex;
    ints = *ints.diagonalize_blocks(eig, num);
    coeff = make_shared<const Matrix>(*tildex * ints);
  }

  tuple<int,int,int,int> nclosed = atommap.num_closed(ga->atoms().front()->name());
  tuple<int,int,int,int> nopen   = atommap.num_open(ga->atoms().front()->name());
  const array<int,4> nclo {{get<0>(nclosed)/2, get<1>(nclosed)/2, get<2>(nclosed)/2, get<3>(nclosed)/2}};
  const array<int,4> nope {{get<0>(nopen),     get<1>(nopen),     get<2>(nopen),     get<3>(nopen)}};
  const int sclosed = (!ga->atoms().front()->use_ecp_basis()) ? get<0>(nclosed)+get<1>(nclosed)+get<2>(nclosed)+get<3>(nclosed)
                    : (get<0>(nclosed)+get<1>(nclosed)+get<2>(nclosed)+get<3>(nclosed) - ga->atoms().front()->ecp_parameters()->ecp_ncore());
  const int sopen   = get<0>(nopen)  +get<1>(nopen)  +get<2>(nopen)  +get<3>(nopen);
  if (sclosed+sopen != ga->nele())
    throw logic_error("Inconsistent nclosed and nopen. See AtomMap");

  auto ocoeff = make_shared<Matrix>(ga->nbasis(), ga->nbasis());
  auto vcoeff = make_shared<Matrix>(ga->nbasis(), ga->nbasis());

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

  auto vden = make_shared<Matrix>(*vcoeff ^ *vcoeff);

  int iter = 0;
  const int maxiter = 100;
  double prev_energy = 0.0;
  for (; iter != maxiter; ++iter) {
    shared_ptr<Matrix> fock = sclosed ? make_shared<Fock<1>>(ga, hcore, nullptr, ocoeff->slice(0,sclosed/2), false/*store*/, true/*rhf*/)
                                      : hcore->copy();
    shared_ptr<Matrix> fock2 = make_shared<Fock<1>>(ga, hcore, vden, ocoeff /*unused*/, false/*store*/, false/*rhf*/, 0.0/*exch*/);
    *fock += *fock2 - *hcore;

    auto aodensity = make_shared<const Matrix>((*ocoeff^*ocoeff)*2.0 + *vden);
    const double energy = ((*hcore+*fock) * *aodensity).trace()*0.5;
    cout << setprecision(10) << energy << endl;

    auto residual = make_shared<const Matrix>(*fock**aodensity**overlap - *overlap**aodensity**fock);
    if (residual->rms() < 1.0e-2 || fabs(energy-prev_energy) < 1.0e-4) break;
    prev_energy = energy;
    fock = diis.extrapolate({fock, residual});

    Matrix ints = *tildex % *fock * *tildex;
    ints = *ints.diagonalize_blocks(eig, num);

    coeff  = make_shared<const Matrix>(*tildex * ints);
    ocoeff = make_shared<Matrix>(ga->nbasis(), sclosed);
    vcoeff = make_shared<Matrix>(ga->nbasis(), ga->nbasis());
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
    vden = make_shared<Matrix>(*vcoeff ^ *vcoeff);
  }
  if (iter == maxiter) throw runtime_error("spin-averaged atomic HF did not converge");

  shared_ptr<Matrix> out = vden;
  if (sclosed) {
    auto c = make_shared<const Coeff>(*ocoeff);
    *out += *c->form_density_rhf(ocoeff->mdim());
  }
  return out;
}
