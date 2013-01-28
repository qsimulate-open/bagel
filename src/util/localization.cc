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

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/util/localization.h>
#include <src/scf/overlap.h>
#include <src/scf/coeff.h>
#include <src/util/jacobi.h>

using namespace bagel;
using namespace std;

/************************************************************************************
* Orthogonalize based on regions - follows the implementation in                    *
*   de Silva, Giebultowski, Korchowiec, PCCP 2011, 14, 546–552                      *
************************************************************************************/
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

  sqrt_S_ = shared_ptr<Matrix>(new Overlap(geom_)); sqrt_S_->sqrt();
  S_inverse_half_ = shared_ptr<Matrix>(new Overlap(geom_)); S_inverse_half_->inverse_half();
}

shared_ptr<Matrix> RegionLocalization::localize_space(shared_ptr<Matrix> density) {
  const int nbasis = geom_->nbasis();

  // Symmetric orthogonalized density matrix
  shared_ptr<Matrix> ortho_density(new Matrix( (*sqrt_S_) * (*density) * (*sqrt_S_) ));

  // transform will hold the eigenvectors of each block
  vector<double> eigenvalues(nbasis, 0.0);
  shared_ptr<Matrix> T = ortho_density->diagonalize_blocks(eigenvalues.data(), sizes_);

  *ortho_density = (*T) % (*ortho_density) * (*T);

  // U matrix will collect the transformations. It should start as the identity.
  shared_ptr<Matrix> U(new Matrix(nbasis, nbasis)); U->add_diag(1.0);

  // Classify each eigenvalue as occupied, mixed, or virtual. All of these classifications may need to be revisited at some point.
  vector<int> occupied, mixed, virt;
  {
    int imo = 0;
    for(auto eig_iter = eigenvalues.begin(); eig_iter != eigenvalues.end(); ++eig_iter, ++imo) {
      if (*eig_iter > 1.5) occupied.push_back(imo);
      else if (*eig_iter > 0.5) mixed.push_back(imo);
      else virt.push_back(imo);
    }
  }

  // Starting rotations
  {
    shared_ptr<JacobiDiag> jacobi(new JacobiDiag(ortho_density,U));

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
  shared_ptr<Matrix> tmp(new Matrix((*S_inverse_half_) * (*T) * (*U)));

  shared_ptr<Matrix> out(new Matrix(nbasis, nbasis));

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

  return out;
}

shared_ptr<const Coeff> RegionLocalization::localize(const int iter, const double thresh) {
  const int nbasis = geom_->nbasis();

  shared_ptr<const Coeff> out;

  shared_ptr<Matrix> closed_coeff = localize_space(coeff_->form_density_rhf(nclosed_));
  if (nact_ == 0) { // In this case, I'm done
    out = shared_ptr<const Coeff>(new const Coeff(*closed_coeff));
  }
  else {
    shared_ptr<Matrix> tmp(new Matrix(nbasis, nbasis));

    vector<double> active_weights(nbasis, 0.0);
    fill_n(active_weights.begin() + nclosed_, nact_, 1.0);
    shared_ptr<Matrix> active_coeff = localize_space(coeff_->form_weighted_density_rhf(nclosed_+nact_, active_weights));

    vector<double> virt_weights(nbasis, 0.0);
    fill(virt_weights.begin() + nclosed_ + nact_, virt_weights.begin() + nbasis, 1.0);
    shared_ptr<Matrix> virt_coeff = localize_space(coeff_->form_weighted_density_rhf(nbasis, virt_weights));

    tmp->copy_block(0, 0, nbasis, nclosed_, closed_coeff);
    tmp->copy_block(0, nclosed_, nbasis, nact_, active_coeff);
    tmp->copy_block(0, nclosed_ + nact_, nbasis, nbasis - (nclosed_ + nact_), virt_coeff);

    out = shared_ptr<const Coeff>(new const Coeff(*tmp));
  }

  return out;
}

void PMLocalization::common_init(shared_ptr<const Geometry> geom) {
  S_ = shared_ptr<Matrix>(new Overlap(geom));

  vector<vector<int>> offsets = geom->offsets();
  for(auto ioffset = offsets.begin(); ioffset != offsets.end(); ++ioffset) {
    const int start = ioffset->front();
    int end;
    if((ioffset+1) == offsets.end()) end = geom->nbasis();
    else end = (ioffset+1)->front();

    atom_bounds_.push_back(make_pair(start, end));
  }
}

shared_ptr<const Coeff> PMLocalization::localize(const int iter, const double thresh) {
  iter_ = iter; thresh_ = thresh;

  shared_ptr<Matrix> new_coeff(new Matrix(*coeff_));
  cout << " === Starting Pipek-Mezey Localization ===" << endl << endl;

  // Localize occupied space
  cout << "  Localizing occupied space" << endl;
  localize_space(new_coeff, 0, nclosed_);

  // If there is an active space, localize it
  if (nact_ != 0) {
    cout << "  Localizing active space" << endl;
    localize_space(new_coeff, nclosed_, nact_);
  }
  else {
    cout << "  No active space to localize" << endl << endl;
  }
    

  // If there is virtual left, localize it
  if ( nvirt_ != 0 ) {
    cout << "  Localizing virtual space" << endl;
    localize_space(new_coeff, nclosed_ + nact_, nvirt_);
  }
  else {
    cout << "  No virtual space to localize" << endl << endl;
  }

  shared_ptr<const Coeff> out(new const Coeff(*new_coeff));
  coeff_ = out;

  return out;
}

void PMLocalization::localize_space(shared_ptr<Matrix> coeff, const int nstart, const int norb) {
  shared_ptr<JacobiPM> jacobi(new JacobiPM(coeff, nstart, norb, S_, atom_bounds_));

  cout << "Iteration          P              dP" << endl;

  double P = calc_P(coeff, nstart, norb);

  cout << setw(3) << 0 << setw(16) << setprecision(10) << P << endl;

  for(int i = 0; i < iter_; ++i) {
    jacobi->sweep();

    double tmp_P = calc_P(coeff, nstart, norb);
    double dP = tmp_P - P;
    cout << setw(3) << i+1 << setw(16) << setprecision(10) << tmp_P
                           << setw(16) << setprecision(10) << dP << endl;
    P = tmp_P;
    if (abs(dP)/P < thresh_) {
      cout << "Converged!" << endl;
      break;
    }
  }
  cout << endl;
}

double PMLocalization::calc_P(shared_ptr<const Matrix> coeff, const int nstart, const int norb) const {
  double out = 0.0;

  const int nbasis = coeff->ndim();
  for(int imo = nstart; imo < (norb + nstart); ++imo) {
    for(auto& ibounds : atom_bounds_) {
      double QA = 0.0;
      for(int iao = ibounds.first; iao < ibounds.second; ++iao) {
        QA += coeff->element(iao,imo) * ddot_(nbasis, coeff->element_ptr(0,imo), 1, S_->element_ptr(0,iao), 1);
      }
      out += QA*QA;
    }
  }

  return out;
}

double PMLocalization::metric() const {
  return calc_P(coeff_, 0, nclosed_);
}
