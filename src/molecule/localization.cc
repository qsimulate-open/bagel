//
// BAGEL - Parallel electron correlation program.
// Filename: localization.h
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

#include <src/molecule/overlap.h>
#include <src/math/jacobi.h>
#include <src/molecule/localization.h>

using namespace bagel;
using namespace std;

OrbitalLocalization::OrbitalLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  const int nclosed, const int nact, const int nvirt) :
  input_(input), geom_(geom), coeff_(coeff->distmatrix()), nclosed_(nclosed), nact_(nact), nvirt_(nvirt) {

  localize_closed_ = input->get<bool>("closed", true);
  localize_active_ = input->get<bool>("active", (nact_ > 0));
  localize_virtual_ = input->get<bool>("virtual", false);

  thresh_ = input->get<double>("thresh", 1.0e-8);
  max_iter_ = input->get<int>("max_iter", 50);
}

/************************************************************************************
* Orthogonalize based on regions - follows the implementation in                    *
*   de Silva, Giebultowski, Korchowiec, PCCP 2011, 14, 546–552                      *
************************************************************************************/
RegionLocalization::RegionLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  vector<int> region_sizes, const int nclosed, const int nact, const int nvirt) :
    OrbitalLocalization(input, geom, coeff, nclosed, nact, nvirt)
{
  common_init(region_sizes);
}

RegionLocalization::RegionLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  const int nclosed, const int nact, const int nvirt) :
    OrbitalLocalization(input, geom, coeff, nclosed, nact, nvirt)
{
  vector<int> region_sizes = input->get_vector<int>("region_sizes");
  if (region_sizes.size() < 2) throw runtime_error("At least two regions should be specified for region localization");
  common_init(region_sizes);
}

void RegionLocalization::common_init(vector<int> sizes) {
  int natom = 0;
  for (auto& isize : sizes) {
    int start = geom_->offset(natom).front();
    int end = (((natom + isize) >= (geom_->natom())) ? geom_->nbasis() : geom_->offset(natom + isize).front() );
    sizes_.push_back(end - start);
    bounds_.push_back(make_pair(start, end));

    natom += isize;
  }

  if (natom != geom_->natom()) throw runtime_error("Improper number of atoms in the defined regions of RegionLocalization");

  shared_ptr<Matrix> lsqS = make_shared<Overlap>(geom_); lsqS->sqrt();
  sqrt_S_ = lsqS->distmatrix();
  shared_ptr<Matrix> lSih = make_shared<Overlap>(geom_); lSih->inverse_half();
  S_inverse_half_ = lSih->distmatrix();
}

shared_ptr<DistMatrix> RegionLocalization::localize_space(shared_ptr<const DistMatrix> density) {
  const int nbasis = geom_->nbasis();
  shared_ptr<const Matrix> den = density->matrix();

  // Symmetric orthogonalized density matrix
  auto ortho_density = make_shared<Matrix>((*sqrt_S_) * (*den) * (*sqrt_S_));

  // transform will hold the eigenvectors of each block
  vector<double> eigenvalues(nbasis, 0.0);
  shared_ptr<Matrix> T = ortho_density->diagonalize_blocks(eigenvalues.data(), sizes_);

  *ortho_density = (*T) % (*ortho_density) * (*T);

  // U matrix will collect the transformations. It should start as the identity.
  auto U = make_shared<Matrix>(nbasis, nbasis); U->unit();

  // Classify each eigenvalue as occupied, mixed, or virtual. All of these classifications may need to be revisited at some point.
  vector<int> region_occupied(sizes_.size(), 0); // Will contain info on number of occ/virt orbitals in each region
  vector<int> occupied, mixed, virt;
  {
    int imo = 0;
    const int nregions = sizes_.size();
    for (int iregion = 0; iregion < nregions; ++iregion) {
      const int isize = sizes_.at(iregion);
      for(int i = 0; i < isize; ++i, ++imo) {
        if (eigenvalues.at(imo) > 1.5) { occupied.push_back(imo); region_occupied.at(iregion) += 1; }
        else if (eigenvalues.at(imo) > 0.5) mixed.push_back(imo);
        else { virt.push_back(imo); }
      }
    }
  }
  region_orbitals_.push_back(region_occupied);

  // Starting rotations
  {
    auto jacobi = make_shared<JacobiDiag>(make_shared<PTree>(), ortho_density->distmatrix(), U->distmatrix());

    for(int& iocc : occupied) {
      for(int& imixed : mixed) jacobi->rotate(iocc, imixed);
      for(int& ivirt : virt) jacobi->rotate(iocc, ivirt);
    }

    for(int& imixed : mixed) {
      for(int& ivirt : virt) jacobi->rotate(imixed, ivirt);
    }
  }

  // TODO: here would be the stage 2 jacobi sweeps. This would be necessary for covalently bound regions.

  // Reorder so that the occupied orbitals come first, separated by fragment
  auto tmp = make_shared<Matrix>((*S_inverse_half_) * (*T) * (*U));
  auto out = make_shared<Matrix>(nbasis, nbasis);
  {
    int imo = 0;
    for(int& iocc : occupied) {
      copy_n(tmp->element_ptr(0,iocc), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
    for(int& imixed : mixed) {
      copy_n(tmp->element_ptr(0,imixed), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
    for(int& ivirt : virt) {
      copy_n(tmp->element_ptr(0,ivirt), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
  }

  return out->distmatrix();
}

shared_ptr<const Matrix> RegionLocalization::localize() {
  const int nbasis = geom_->nbasis();

  shared_ptr<const Matrix> closed = localize_space(coeff_->form_density_rhf(nclosed_))->matrix();

  if (!localize_active_) {
    coeff_ = closed->distmatrix();
    return closed;
  }
  else {
    auto out = make_shared<Matrix>(*coeff_);
    copy_n(closed->element_ptr(0, 0), nclosed_ * nbasis, out->element_ptr(0, 0));

    shared_ptr<Matrix> active = localize_space(coeff_->form_density_rhf(nact_, nclosed_))->matrix();
    copy_n(active->element_ptr(0, 0), nact_ * nbasis, out->element_ptr(0, nclosed_));

    if (localize_virtual_) {
      shared_ptr<Matrix> virt = localize_space(coeff_->form_density_rhf(nvirt_, nclosed_ + nact_))->matrix();
      copy_n(virt->element_ptr(0, 0), nvirt_ * nbasis, out->element_ptr(0, nclosed_ + nact_));
    }

    return out;
  }
}

/************************************************************************************
* Pipek-Mezey Localization                                                          *
************************************************************************************/
PMLocalization::PMLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  const int nclosed, const int nact, const int nvirt) : OrbitalLocalization(input, geom, coeff, nclosed, nact, nvirt)
{
  common_init();
}

void PMLocalization::common_init() {
  shared_ptr<Matrix> lS = make_shared<Overlap>(geom_);
  S_ = lS->distmatrix();

  vector<vector<int>> offsets = geom_->offsets();
  int nbasis = 0;
  for(auto& atom : geom_->atoms()) {
    const int start = nbasis;
    nbasis +=  atom->nbasis();
    if (start != nbasis)
      atom_bounds_.emplace_back(start, nbasis);
  }
  assert(nbasis == geom_->nbasis());
}

shared_ptr<const Matrix> PMLocalization::localize() {
  auto out = make_shared<DistMatrix>(*coeff_);

  cout << " === Starting Pipek-Mezey Localization ===" << endl << endl;

  // Localize closed space
  if (localize_closed_) {
    cout << "  Localizing closed (occupied) space" << endl;
    localize_space(out, 0, nclosed_);
  }
  else {
    cout << "  Skipping closed (occupied) space" << endl << endl;
  }

  // If there is an active space, localize it
  if (localize_active_) {
    cout << "  Localizing active space" << endl;
    localize_space(out, nclosed_, nact_);
  }
  else {
    cout << "  Skipping active space" << endl << endl;
  }

  // If there is virtual left, localize it
  if ( localize_virtual_ ) {
    cout << "  Localizing virtual space" << endl;
    localize_space(out, nclosed_ + nact_, nvirt_);
  }
  else {
    cout << "  No virtual space to localize" << endl << endl;
  }

  coeff_ = out;

  return out->matrix();
}

void PMLocalization::localize_space(shared_ptr<DistMatrix> coeff, const int nstart, const int norb) {
  auto jacobi = make_shared<JacobiPM>(input_, coeff, nstart, norb, S_, atom_bounds_);

  cout << "Iteration          P              dP" << endl;

  double P = calc_P(coeff, nstart, norb);

  cout << setw(3) << 0 << setw(16) << setprecision(10) << P << endl;

  for(int i = 0; i < max_iter_; ++i) {
    jacobi->sweep();

    double tmp_P = calc_P(coeff, nstart, norb);
    double dP = tmp_P - P;
    cout << setw(3) << i+1 << setw(16) << setprecision(10) << tmp_P
                           << setw(16) << setprecision(10) << dP << endl;
    P = tmp_P;
    if (fabs(dP) < thresh_) {
      cout << "Converged!" << endl;
      break;
    }
  }
  cout << endl;
}

double PMLocalization::calc_P(shared_ptr<const DistMatrix> coeff, const int nstart, const int norb) const {
  const int nbasis = coeff->ndim();

  double out = 0.0;
  auto mos = make_shared<DistMatrix>(nbasis, norb);

#ifdef HAVE_SCALAPACK
  pdgemm_("N", "N", nbasis, norb, nbasis, 1.0, S_->local().get(), 1, 1, S_->desc().get(),
                                                 coeff->local().get(), 1, nstart + 1, coeff->desc().get(),
                                            0.0, mos->local().get(), 1, 1, mos->desc().get());
#else
  dgemm_("N", "N", nbasis, norb, nbasis, 1.0, S_->data(), nbasis, coeff->element_ptr(0, nstart), nbasis, 0.0, mos->data(), nbasis);
#endif

  const int natom = atom_bounds_.size();

  auto P_A = make_shared<DistMatrix>(norb, norb);

  for (auto& ibounds : atom_bounds_) {
    const int natombasis = ibounds.second - ibounds.first;

#ifdef HAVE_SCALAPACK
    pdgemm_("T", "N", norb, norb, natombasis, 1.0, mos->local().get(), ibounds.first + 1, 1, mos->desc().get(),
                                                         coeff->local().get(), ibounds.first + 1, nstart + 1, coeff->desc().get(),
                                                    0.0, P_A->local().get(), 1, 1, P_A->desc().get());
#else
    dgemm_("T", "N", norb, norb, natombasis, 1.0, mos->element_ptr(ibounds.first, 0), nbasis,
                              coeff->element_ptr(ibounds.first, nstart), nbasis, 0.0, P_A->data(), norb);
#endif

    shared_ptr<const Matrix> localPA = P_A->matrix();

    for (int imo = 0; imo < norb; ++imo)
      out += localPA->element(imo, imo) * localPA->element(imo, imo);
  }

  return out;
}

double PMLocalization::metric() const {
  return calc_P(coeff_, 0, nclosed_);
}
