//
// BAGEL - Parallel electron correlation program.
// Filename: ericompute.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

  void set_data (const double* pos, const double exp, const int* ang_mom, const std::vector<double> field);
  void change_angular (const int ax, const int ay, const int az);
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


class nucleus {

  public:
  int charge;
  double exponent;
  double position[3];

  //constructor
  nucleus (int Z, std::vector<double> coords, double exp) {
    if (Z<1) throw std::runtime_error ("Nucleus should have a positive charge.  Check that inputs are being read properly.");
    charge = Z;
    for (int i=0; i!=3; i++) position[i] = coords[i];
    exponent = exp;
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
ryan::polynomial<std::complex<double>> get_ERI_Ix (const int dimension, const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
ryan::polynomial<std::complex<double>> get_ERI_III (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);


// Setting up MOs, etc. before computing integrals
std::pair<std::vector<atomic_orbital>,std::vector<molecular_orbital>> prepare_orbitals (int nbasis_contracted, bool normalize_basis, bool scale_input, bool orthogonalize,
    std::vector<double> field, std::vector<double> positions, std::vector<int> angular, std::vector<double> exponents,
    std::vector<double> contraction_coefficients, std::vector<int> nprimitive, std::vector<std::complex<double>> orbital1,
    std::vector<std::complex<double>> orbital2, std::vector<std::complex<double>> orbital3, std::vector<std::complex<double>> orbital4);


// Used in computation of ERI
std::pair<std::complex<double>,std::complex<double>> compute_eri_ssss(const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
std::complex<double> get_eri_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_);
std::complex<double> compute_eri (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field);

// Used in computation of NAI
std::pair<std::complex<double>,std::complex<double>> compute_ss(const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_);
std::complex<double> get_nai_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_);
std::complex<double> compute_nai (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, std::vector<nucleus> nuclei);
std::complex<double> compute_finite_nai (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, std::vector<nucleus> nuclei);
std::complex<double> get_finite_nai_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_);

ryan::polynomial<std::complex<double>> get_NAI_III (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_);
ryan::polynomial<std::complex<double>> get_NAI_Ix (const int dimension, const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_);

// Used in computation of kinetic energy
std::complex<double> kinetic_MO (std::vector<double> field, molecular_orbital A_, molecular_orbital B_, std::vector<atomic_orbital> basis);
std::complex<double> kinetic (std::vector<double> field, atomic_orbital A_, atomic_orbital B_);
std::vector<std::complex<double>> momentum_MO (std::vector<double> field, molecular_orbital A_, molecular_orbital B_, std::vector<atomic_orbital> basis);
std::vector<std::complex<double>> momentum (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_);

// Used for small component integrals
std::complex<double> compute_smallnai (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, std::vector<nucleus> nuclei, const int ia, const int ib);
std::complex<double> compute_small_finitenai (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, std::vector<nucleus> nuclei, const int ia, const int ib);
std::complex<double> get_smallnai_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_, const int ia, const int ib);
std::complex<double> get_small_finitenai_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, nucleus C_, const int ia, const int ib);
std::complex<double> compute_smalloverlap (std::vector<double> field, molecular_orbital A_, molecular_orbital B_, std::vector<atomic_orbital> basis, const int ia, const int ib);
std::complex<double> get_smalloverlap_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, const int ia, const int ib);
std::complex<double> compute_smalleri (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, const int ia, const int ib);
std::complex<double> get_smalleri_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_, const int ia, const int ib);
std::complex<double> compute_mixederi (std::vector<atomic_orbital> basis, std::vector<molecular_orbital> input, std::vector<double> field, const int ia, const int ib);
std::complex<double> get_mixederi_matrix_element (const std::vector<double> field, atomic_orbital A_, atomic_orbital B_, atomic_orbital C_, atomic_orbital D_, const int ia, const int ib);

}

#endif
