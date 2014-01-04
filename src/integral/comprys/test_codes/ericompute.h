//
// BAGEL - Parallel electron correlation program.
// Filename: ericompute.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan Reynolds <rreynoldschem@u.northwestern.edu>
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

// The functions defined here allow us to compute the ERI (and overlap integral) for London orbitals in a straightforward but inefficient way.
// Their main purpose is to aid in debugging of Bagel's actual integral codes.

#ifndef __ERICOMPUTE_H
#define __ERICOMPUTE_H

#include <array>
#include <vector>
#include <cmath>
#include <complex>
#include "../../../molecule/shell.h"
#include "polynomial.h"

namespace ryan {

constexpr double pi = 3.1415926535897932385;

class atomic_orbital {

  public:
  double prefactor;
  double exponent;
  double position[3];
  double vector_potential[3];
  int angular_momentum[3];

  atomic_orbital () {
    prefactor = 0.0;

  }
  void set_data (const double* pos, const double exp, const int* ang_mom, bool norm, const std::vector<double> field);
};


class molecular_orbital {

  public:
  std::vector<std::complex<double>> coefficient;

  void reassign (std::vector<std::complex<double>> coeffs) {
    coefficient = coeffs;
  }

  //constructor
  molecular_orbital (std::vector<std::complex<double>> coeffs) {
    coefficient = coeffs;
  }
};


// Overlap integral for atomic London orbitals
std::complex<double> overlap_Ix (const int dimension, const std::vector<double> field, atomic_orbital A_, atomic_orbital B_);
std::complex<double> overlap (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_);


// Needed for orthogonalization of basis set
std::complex<double> overlap_MO (std::vector<double> field, molecular_orbital A_, molecular_orbital B_, std::vector<atomic_orbital> basis);
molecular_orbital projection_MO (std::vector<double> field, molecular_orbital A_, molecular_orbital B_, std::vector<atomic_orbital> basis);
molecular_orbital add_MOs (molecular_orbital A_, molecular_orbital B_);
molecular_orbital scalar_MO (std::complex<double> scalar, molecular_orbital A_);
molecular_orbital subtract_MOs (molecular_orbital A_, molecular_orbital B_);
std::vector<molecular_orbital> orthogonalize_basis (std::vector<double> field, std::vector<molecular_orbital> old_basis, std::vector<atomic_orbital> basis);


// ERI Recurrence Relations
ryan::polynomial<std::complex<double>> get_Ix (const int dimension, const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
ryan::polynomial<std::complex<double>> get_III (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);


// Used in computation of ERI
std::pair<std::complex<double>,std::complex<double>> compute_ssss(const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
std::complex<double> get_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
std::complex<double> compute_eri (int nbasis_contracted, bool normalize_basis, bool scale_input, bool orthogonalize, std::vector<double> field,
    std::vector<double> positions, std::vector<int> angular, std::vector<double> exponents, std::vector<double> contraction_coefficients, std::vector<int> nprimitive,
    std::vector<std::complex<double>> orbital1, std::vector<std::complex<double>> orbital2,
    std::vector<std::complex<double>> orbital3, std::vector<std::complex<double>> orbital4);

std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_ERI (const std::array<std::shared_ptr<const bagel::Shell>,4>& basisinfo);

}

#endif
