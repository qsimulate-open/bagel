//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bagel_interface.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// The functions declared here help to connect the inefficent comparison function to Bagel's actual integral evaluation codes.

#ifndef __BAGEL_INTERFACE_H
#define __BAGEL_INTERFACE_H

#include <vector>
#include <complex>
#include "src/molecule/shell.h"
#include "src/molecule/molecule.h"

namespace test {

std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_ERI (const std::array<std::shared_ptr<const bagel::Shell>,4>& basisinfo, const std::array<double,3> magnetic_field);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_NAI (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::shared_ptr<const bagel::Molecule> mol, const std::array<double,3> magnetic_field);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_finiteNAI(const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::shared_ptr<const bagel::Atom> atom, const std::array<double,3> magnetic_field);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_overlap (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::array<double,3> magnetic_field);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_kinetic (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::array<double,3> magnetic_field);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_momentum (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::array<double,3> magnetic_field);

std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_smallNAI (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::shared_ptr<const bagel::Molecule> mol, const std::array<double,3> magnetic_field, const int ia, const int ib);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_small_finiteNAI (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::shared_ptr<const bagel::Molecule> mol, const std::array<double,3> magnetic_field, const int ia, const int ib);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_smalloverlap (const std::array<std::shared_ptr<const bagel::Shell>,2>& basisinfo, const std::array<double,3> magnetic_field, const int ia, const int ib);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_smallERI (const std::array<std::shared_ptr<const bagel::Shell>,4>& basisinfo, const std::array<double,3> magnetic_field, const int ia, const int ib);
std::vector<std::pair<std::vector<int>,std::complex<double>>> get_comparison_mixedERI (const std::array<std::shared_ptr<const bagel::Shell>,4>& basisinfo, const std::array<double,3> magnetic_field, const int ia, const int ib);

}

#endif
