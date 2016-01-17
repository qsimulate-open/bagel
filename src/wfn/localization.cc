//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: localization.h
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

#include <algorithm>
#include <src/mat1e/overlap.h>
#include <src/util/math/jacobi.h>
#include <src/wfn/localization.h>

using namespace bagel;
using namespace std;

OrbitalLocalization::OrbitalLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  vector<pair<int,int>> subspaces) :
  input_(input), geom_(geom), coeff_(coeff), orbital_subspaces_(subspaces)
{
  // extra subspaces defined as a list of 2-member lists: [ [first, last], [first, last], ... ]
  //   the entire range will be localized (from first to last, including last)
  shared_ptr<const PTree> subs = input->get_child_optional("subspaces");
  if (subs) {
    for (auto& sp : *subs) {
      array<int,2> tmp = sp->get_array<int,2>("");
      // accept orbital numbers with 1-based numbering
      orbital_subspaces_.emplace_back(tmp[0]-1, tmp[1]);
    }
  }
}

OrbitalLocalization::OrbitalLocalization(shared_ptr<const PTree> input, shared_ptr<const Reference> ref) :
  OrbitalLocalization(input, ref->geom(), ref->coeff(), vector<pair<int,int>>()) {
  if (input->get<bool>("occupied",true))
    orbital_subspaces_.emplace_back(0,ref->nclosed());
  if (ref->nact() != 0 && input->get<bool>("active",true))
    orbital_subspaces_.emplace_back(ref->nclosed(), ref->nclosed()+ref->nact());
  if (input->get<bool>("virtual", false))
    orbital_subspaces_.emplace_back(ref->nclosed()+ref->nact(),ref->nclosed()+ref->nact()+ref->nvirt());

  // if there, take eigenvalues out of Reference to reorder subspaces later on
  // TODO: would more flexibility in defining the reordering criterion be helpful?
  if (!ref->eig().empty()) {
    diagonals_ = ref->eig();
  }
}

shared_ptr<Matrix> OrbitalLocalization::localize() {
  auto out = make_shared<Matrix>(*coeff_);
  set<int> localized;

  shared_ptr<Matrix> ordering;
  if (!diagonals_.empty()) {
    ordering = out->copy();
    for (int i = 0; i < ordering->mdim(); ++i) {
      const double e = std::sqrt(diagonals_(i));
      for_each(ordering->element_ptr(0,i), ordering->element_ptr(0,i+1), [&e](double& a){ return a *= e; });
    }
    ordering = make_shared<Matrix>(*ordering ^ *ordering);
  }

  // localize within each subspace
  for (auto& sp : orbital_subspaces_) {
    const int start = sp.first;
    const int fence = sp.second;

    cout << "  *** Localizing orbitals " << start+1 << " -- " << fence << " ***" << endl;

    for (int iter = start; iter < fence; ++iter)
      if (!localized.insert(iter).second) cout << "WARNING: Orbital " << iter << " has already been localized." << endl;

    shared_ptr<Matrix> loc_subspace = localize_space(out->slice_copy(start, fence));
    if (ordering) {
      Matrix tmp = *loc_subspace % *ordering * *loc_subspace;
      multimap<double, int> to_copy;
      for (int i = 0; i < tmp.ndim(); ++i) to_copy.emplace(tmp(i,i), i);
      for (auto& i : to_copy) copy_n(loc_subspace->element_ptr(0,i.second), loc_subspace->ndim(), out->element_ptr(0,i.second+start));
    }
    else {
      out->copy_block(0, start, out->ndim(), (fence-start), loc_subspace);
    }
  }

  // is this necessary for anything but the test suite?
  coeff_ = out;
  return out;
}

/************************************************************************************
* Orthogonalize based on regions - follows the implementation in                    *
*   de Silva, Giebultowski, Korchowiec, PCCP 2011, 14, 546–552                      *
************************************************************************************/
RegionLocalization::RegionLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  vector<pair<int, int>> subspaces, vector<int> region_sizes) :
    OrbitalLocalization(input, geom, coeff, subspaces)
{
  common_init(region_sizes);
}

RegionLocalization::RegionLocalization(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, vector<int> region_sizes) :
  OrbitalLocalization(input, ref)
{
  common_init(region_sizes);
}

void RegionLocalization::common_init(vector<int> sizes) {
  cout << " ======    Region Localization    ======" << endl;
  vector<int> input_sizes = input_->get_vector<int>("region_sizes");
  sizes.insert(sizes.end(), input_sizes.begin(), input_sizes.end());
  if (sizes.size() < 2) throw runtime_error("At least two regions should be specified for region localization");

  int natom = 0;
  for (auto& isize : sizes) {
    int start = geom_->offset(natom).front();
    int end = (((natom + isize) >= (geom_->natom())) ? geom_->nbasis() : geom_->offset(natom + isize).front() );
    sizes_.push_back(end - start);
    region_bounds_.emplace_back(start, end);

    natom += isize;
  }

  if (natom != geom_->natom()) throw runtime_error("Improper number of atoms in the defined regions of RegionLocalization");

  auto S = make_shared<Overlap>(geom_);
  sqrt_S_ = S->copy(); sqrt_S_->sqrt();
  S_inverse_half_ = S->copy();
  if (!S_inverse_half_->inverse_half()) {
    throw runtime_error("Region Localization does not handle linear dependencies. Use PM localization.");
  }
}

shared_ptr<Matrix> RegionLocalization::localize_space(shared_ptr<const Matrix> coeffs) {
  Matrix den = *coeffs ^ *coeffs;
  den.scale(2.0);

  const int nbasis = geom_->nbasis();

  // Symmetric orthogonalized density matrix
  auto ortho_density = make_shared<Matrix>(*sqrt_S_ % den * *sqrt_S_);

  // transform will hold the eigenvectors of each block
  VectorB eigenvalues(nbasis);
  Matrix T = *ortho_density->diagonalize_blocks(eigenvalues, sizes_);

  *ortho_density = T % *ortho_density * T;

  // U matrix will collect the transformations. It should start as the identity.
  auto U = make_shared<Matrix>(nbasis, nbasis); U->unit();

  // Classify each eigenvalue as occupied, mixed, or virtual. All of these classifications may need to be revisited at some point.
  vector<int> occupied, mixed, virt;
  {
    int imo = 0;
    const int nregions = sizes_.size();
    for (int iregion = 0; iregion < nregions; ++iregion) {
      const int isize = sizes_[iregion];
      for(int i = 0; i < isize; ++i, ++imo) {
        if (eigenvalues[imo] > 1.5)
          occupied.push_back(imo);
        else if (eigenvalues[imo] > 0.5)
          mixed.push_back(imo);
        else
          virt.push_back(imo);
      }
    }
  }

  // Starting rotations
  {
    auto jacobi = make_shared<JacobiDiag>(make_shared<PTree>(), ortho_density, U);

    for(int& iocc : occupied) {
      for(int& imixed : mixed) jacobi->rotate(iocc, imixed);
      for(int& ivirt : virt) jacobi->rotate(iocc, ivirt);
    }

    for(int& imixed : mixed) {
      for(int& ivirt : virt) jacobi->rotate(imixed, ivirt);
    }

    if ( !mixed.empty() )
      cout << "WARNING! Localization between bound regions not well tested." << endl;
    for (auto i = mixed.begin(); i != mixed.end(); ++i) {
      for (auto j = mixed.begin(); j != i; ++j)
        jacobi->rotate(*j,*i);
    }
  }

  // Reorder so that the occupied orbitals come first, separated by fragment
  Matrix localized_coeffs = *S_inverse_half_ * T * *U;
  auto out = make_shared<Matrix>(nbasis, occupied.size() + mixed.size()/2);
  {
    int imo = 0;
    for(int& iocc : occupied) {
      copy_n(localized_coeffs.element_ptr(0,iocc), nbasis, out->element_ptr(0,imo));
      ++imo;
    }
    for(int& imixed : mixed) {
      if (ortho_density->element(imixed,imixed) > 1.5) {
        copy_n(localized_coeffs.element_ptr(0,imixed), nbasis, out->element_ptr(0,imo));
        ++imo;
      }
    }
    if (imo != out->mdim())
      throw runtime_error("Unexpected number of orbitals in region localization");
  }

  return out;
}

/************************************************************************************
* Pipek-Mezey Localization                                                          *
************************************************************************************/
PMLocalization::PMLocalization(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom, shared_ptr<const Matrix> coeff,
  vector<pair<int, int>> subspaces, vector<int> sizes) : OrbitalLocalization(input, geom, coeff, subspaces)
{
  common_init(sizes);
}

PMLocalization::PMLocalization(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, vector<int> sizes)
  : OrbitalLocalization(input, ref)
{
  common_init(sizes);
}

void PMLocalization::common_init(vector<int> sizes) {
  cout << " ======    Pipek-Mezey Localization    ======" << endl;

  max_iter_ = input_->get<int>("max_iter", 50);
  thresh_ = input_->get<double>("thresh", 1.0e-6);
  lowdin_ = input_->get<bool>("lowdin", true);

  cout << endl << "  Localization threshold: " << setprecision(2) << setw(6) << scientific << thresh_ << endl << endl;

  S_ = make_shared<Overlap>(geom_);
  if (lowdin_) S_->sqrt();

  string localization_type = input_->get<string>("type", "atomic");

  int nbasis = 0;

  if (localization_type == "atomic") {
    for(auto& atom : geom_->atoms()) {
      const int start = nbasis;
      nbasis +=  atom->nbasis();
      if (start != nbasis)
        region_bounds_.emplace_back(start, nbasis);
    }
  }
  else if (localization_type == "region") {
    if (input_->get_child_optional("region_sizes")) {
      vector<int> input_sizes = input_->get_vector<int>("region_sizes");
      sizes.insert(sizes.end(), input_sizes.begin(), input_sizes.end());
    }
    int natoms = 0;
    for (int& region : sizes) {
      const int atomstart = natoms;
      const int basisstart = nbasis;
      for (int atom = atomstart; atom < atomstart + region; ++atom)
        nbasis += geom_->atoms()[atom]->nbasis();

      natoms += region;
      if (basisstart != nbasis)
        region_bounds_.emplace_back(basisstart, nbasis);
    }
    if ( natoms != count_if(geom_->atoms().begin(), geom_->atoms().end(), [](const shared_ptr<const Atom> a){return !a->dummy();}) ) {
      throw logic_error("All atoms must be assigned to regions");
    }
  }
  else {
    throw logic_error("Unrecognized PM localization type");
  }
  assert(nbasis == geom_->nbasis());
}

shared_ptr<Matrix> PMLocalization::localize_space(shared_ptr<const Matrix> coeff) {
  Timer pmtime;
  auto out = make_shared<Matrix>(*coeff);
  const int norb = out->mdim();

  auto jacobi = make_shared<JacobiPM>(input_, out, 0, norb, S_, region_bounds_, lowdin_);

  //cout << "iteration            P_A^2                delta P_A^2" << endl;
  cout << setw(6) << "iter" << setw(20) << "P_A^2" << setw(27) << "delta P_A^2" << setw(22) << "time" << endl;
  cout << "----------------------------------------------------------------------------------------------" << endl;

  double P = calc_P(out, 0, norb);

  cout << setw(5) << 0 << fixed << setw(24) << setprecision(10) << P << endl;

  for(int i = 0; i < max_iter_; ++i) {
    jacobi->sweep();

    double tmp_P = calc_P(out, 0, norb);
    double dP = tmp_P - P;
    cout << setw(5) << i+1 << fixed << setw(24) << setprecision(10) << tmp_P
                           << fixed << setw(24) << setprecision(10) << dP
                           << fixed << setw(24) << setprecision(6)  << pmtime.tick() << endl;
    P = tmp_P;
    if (fabs(dP) < thresh_) {
      cout << "Converged!" << endl;
      break;
    }
  }
  cout << endl;

  return out;
}

double PMLocalization::calc_P(shared_ptr<const Matrix> coeff, const int nstart, const int norb) const {
  const int nbasis = coeff->ndim();

  double out = 0.0;
  auto mos = make_shared<Matrix>(nbasis, norb);

  dgemm_("N", "N", nbasis, norb, nbasis, 1.0, S_->data(), nbasis, coeff->element_ptr(0, nstart), nbasis, 0.0, mos->data(), nbasis);

  auto P_A = make_shared<Matrix>(norb, norb);

  for (auto& ibounds : region_bounds_) {
    const int natombasis = ibounds.second - ibounds.first;

    if (lowdin_)
      dgemm_("T", "N", norb, norb, natombasis, 1.0, mos->element_ptr(ibounds.first, 0), nbasis,
                              mos->element_ptr(ibounds.first, 0), nbasis, 0.0, P_A->data(), norb);
    else
      dgemm_("T", "N", norb, norb, natombasis, 1.0, mos->element_ptr(ibounds.first, 0), nbasis,
                              coeff->element_ptr(ibounds.first, nstart), nbasis, 0.0, P_A->data(), norb);

    for (int imo = 0; imo < norb; ++imo)
      out += P_A->element(imo, imo) * P_A->element(imo, imo);
  }

  return std::sqrt(out/static_cast<double>(norb));
}

double PMLocalization::metric() const {
  return calc_P(coeff_, 0, geom_->nele()/2);
}
