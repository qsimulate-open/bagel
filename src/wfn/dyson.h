//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dyson.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Alexander Humeniuk
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


#ifndef __BAGEL_WFN_DYSON_H
#define __BAGEL_WFN_DYSON_H

#include <set>
#include <bitset>

#include <src/wfn/reference.h>
#include <src/util/math/vectorb.h>

namespace bagel {

  const std::string indent = "  ";
  
class DysonOrbitals {
  protected:
    // Hang on to the input for convenience
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> geom_;
    // initial wavefunctions (with N electrons)
    std::shared_ptr<const Reference> refI_;
    // final ionized wavefunctions (with N-1 electrons)
    std::shared_ptr<const Reference> refF_;
    // indexes of initial states
    std::vector<int> initial_states_;
    // indexes of final states
    std::vector<int> final_states_;
    // threshold for neglecting products of CI coefficients
    double thresh_;
    // Dyson orbitals are written to this molden file.
    std::string molden_file_;
    
    // AO overlap matrix
    Overlap Sao_;
    // overlap matrix between alpha and beta molecular orbitals
    Matrix SmoA_;
    Matrix SmoB_;

    // number of AOs
    int nao_;
    // number of MOs
    int nmo_;
    // number of ionization channels = (nr. initial states) x (nr. final states)
    int nchan_;
    // MO coefficients of Dyson orbitals (matrix with dimension  nao x nchan_)
    Matrix coeff_;
    // norms of Dyson orbitals
    VectorB norms_;
    // ionization energies 
    VectorB energies_;
    
  private:
    bool consistency_checks();
    void mo_overlaps();
    void ci_dyson();
    VectorB minors(std::shared_ptr<const Matrix> mat);
    VectorB slater_dyson(std::bitset<nbit__> bita1, std::bitset<nbit__> bitb1,
			 std::shared_ptr<const Determinants> det1, int nclosed1,
			 std::bitset<nbit__> bita2, std::bitset<nbit__> bitb2,
			 std::shared_ptr<const Determinants> det2, int nclosed2);

  public:
    DysonOrbitals(std::shared_ptr<const PTree> input);
    void compute();
    void print_results();
    void write_molden();
    // access functions
    int nchan() const { return nchan_ };
    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
    std::shared_ptr<const VectorB> norms() const { return norms_; }
    std::shared_ptr<const VectorB> energies() const { return energies_; }
};

}

#endif
