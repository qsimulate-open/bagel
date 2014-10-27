//
// BAGEL - Parallel electron correlation program.
// Filename: pdfdist.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#ifndef __SRC_PERIODIC_PDFDIST_H
#define __SRC_PERIODIC_PDFDIST_H

#include <src/periodic/pdfdist_ints.h>
#include <src/periodic/pdata.h>

namespace bagel {

/// basically an array of DFDist for periodic systems
class PDFDist {
  protected:
    /// lattice vectors in direct space
    std::vector<std::array<double, 3>> lattice_vectors_;

    const size_t nbasis_;
    const size_t naux_;

    bool serial_;

    ///  dfdist contains 3-index integrals (r0 sL'|aL) sum over L
    std::vector<std::shared_ptr<PDFDist_ints>> dfdist_;

    /// 2-index integrals (i|j_L)^{-1} (sum over L)
    std::shared_ptr<Matrix> data2_;
    void pcompute_2index(const std::vector<std::shared_ptr<const Shell>>& ashell, const double throverlap);

    /// normalised 1-index auxiliary charge <i|.>
    std::shared_ptr<VectorB> data1_;
    void compute_aux_charge(const std::vector<std::shared_ptr<const Shell>>& ashell);

    /// projection matrix P_{ij} = <i|.><.|j>
    std::shared_ptr<const Matrix> projector_;

  public:
    PDFDist(std::vector<std::array<double, 3>> lattice_vectors, const int nbasis, const int nauxbasis,
            const std::vector<std::shared_ptr<const Atom>>& atoms_cell0,
            const std::vector<std::shared_ptr<const Atom>>& aux_atoms,
            const double thresh, const bool serial = false, const std::shared_ptr<Matrix> data2 = nullptr);

    std::vector<std::array<double, 3>> lattice_vectors() const { return lattice_vectors_; }
    int ncell() const { return lattice_vectors_.size(); }

    std::vector<std::shared_ptr<PDFDist_ints>> dfdist() const { return dfdist_; }
    std::shared_ptr<PDFDist_ints> dfdist(const int L) const { return dfdist_[L]; }

    std::shared_ptr<const Matrix> data2() const { return data2_; }
    std::shared_ptr<const VectorB> data1() const { return data1_; }

    size_t nbasis() const { return nbasis_; }
    size_t naux() const { return naux_; }
    bool serial() const { return serial_; }

    /// compute J_{rs} operator (r0 and sL) for all L, given density matrices D_{rs}^L in AO basis
    std::shared_ptr<VectorB> pcompute_coeff(const std::shared_ptr<const PData> density) const;
    std::shared_ptr<PData>   pcompute_Jop_from_coeff(std::shared_ptr<const VectorB> coeff) const;
    std::shared_ptr<PData>   pcompute_Jop(const std::shared_ptr<const PData> density) const;
};

}

#endif
