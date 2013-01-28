//
// BAGEL - Parallel electron correlation program.
// Filename: pmp2.h
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

#ifndef __src_pmp2_pmp2_h
#define __src_pmp2_pmp2_h

#include <stddef.h>
#include <vector>
#include <tuple>
#include <memory>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pgeometry.h>
#include <src/util/pfile.h>
#include <src/util/pcompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/rysint/eribatch.h>
#include <src/rysint/slaterbatch.h>

namespace bagel {

class PMP2 {
  typedef std::shared_ptr<PMatrix1e> RefMatrix;
  typedef std::shared_ptr<PCoeff> RefCoeff;
  typedef std::shared_ptr<PMOFile<std::complex<double>>> RefMOFile;

  protected:
    // Shared pointer for geometry.
    const std::shared_ptr<PGeometry> geom_;
    // PGeomery that has a union of OBS and CABS (in this order)
    // Will be initialized in generate_cabs()
    std::shared_ptr<PGeometry> union_geom_;

    // Coefficients for MO.
    RefCoeff coeff_;

    // Coefficients for CABS (OBS / auxiliary part respectively)
    RefCoeff cabs_obs_;
    RefCoeff cabs_aux_;
    RefCoeff coeff_cabs_;
    RefMatrix coeff_entire_;
    RefMatrix ri_entire_;

    // HF orbital energies
    const std::vector<double> eig_;

    // AO integrals used in several places
    std::shared_ptr<PCompFile<ERIBatch>> eri_obs_;
    std::shared_ptr<PCompFile<SlaterBatch>> stg_;
    std::shared_ptr<PCompFile<SlaterBatch>> stg2_;
    std::shared_ptr<PCompFile<SlaterBatch>> yp_;
    std::shared_ptr<PCompCABSFile<ERIBatch>> eri_cabs_;
    std::shared_ptr<PCompCABSFile<SlaterBatch>> stg_cabs_;
    std::shared_ptr<PCompCABSFile<SlaterBatch>> stg_cabs2_;
    RefMOFile eri_ii_pp_;
    RefMOFile eri_ii_Ai_;
    RefMOFile stg_ii_pp_;
    RefMOFile stg_ii_Ai_;
    RefMOFile yp_ii_ii_;
    RefMOFile stg2_ii_ii_;

    // Target tensors
    RefMOFile X_;
    RefMOFile V_;
    RefMOFile B_;

    // Some misc constants
    bool use_hy2_;
    int nbasis_;
    int nfrc_;
    int nocc_act_;
    int nocc_;
    int nvir_;
    int ncabs_;
    size_t noovv_;

    std::pair<RefCoeff, RefCoeff> generate_CABS();
    std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_RI();
    std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_fock_weighted_RI() const;

    // Hartree and exchange matrix
    RefMatrix ao_hJ_;
    RefMatrix ao_K_;
    RefMatrix hJ_obs_obs_;
    RefMatrix hJ_obs_cabs_;
    RefMatrix hJ_cabs_obs_;
    RefMatrix hJ_cabs_cabs_;
    RefMatrix K_obs_obs_;
    RefMatrix K_obs_cabs_;
    RefMatrix K_cabs_obs_;
    RefMatrix K_cabs_cabs_;
    RefMatrix fock_obs_obs_;
    RefMatrix fock_obs_cabs_;
    RefMatrix fock_cabs_obs_;
    RefMatrix fock_cabs_cabs_;
    // Hartree builder
    const std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_hJ();
    // Exchange builder
    const std::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> generate_K();

    RefMatrix coulomb_runtime_OBS() const;
    RefMatrix exchange_runtime_OBS() const;
    RefMatrix coulomb_runtime(const bool) const;
    RefMatrix exchange_runtime(const bool) const;

    void fill_in_cabs_matices();

  public:
    PMP2(const std::shared_ptr<PGeometry>, const std::shared_ptr<PCoeff>,
         const std::vector<double>, std::shared_ptr<PCompFile<ERIBatch>>, const bool);
    ~PMP2();

    void compute();
    void compute_conv_mp2();

};

}

#endif

