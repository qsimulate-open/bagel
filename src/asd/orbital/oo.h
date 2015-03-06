//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/oo.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#ifndef __BAGEL_ASD_OO_H
#define __BAGEL_ASD_OO_H

#include <src/wfn/reference.h>
#include <src/asd/dimer/dimer.h>
#include <src/asd/orbital/rotfile.h>
#include <src/df/dfblock.h>

namespace bagel {

class ASD_OO : public Method, public std::enable_shared_from_this<ASD_OO> {

  protected:
    // some internal information
    int nocc_; // sum of nact_ + nclosed_
    int nclosed_;
    int nact_;
    int nvirt_;
    int nbasis_; // number of MO orbitals
    int nstate_;
    int max_iter_;
    double thresh_;
    double precond_;

    int nactA_;
    int nactB_;
    std::array<int,3> rasA_;
    std::array<int,3> rasB_;

    std::shared_ptr<Dimer> dimer_;
    std::shared_ptr<PTree> asdinput_;

    std::shared_ptr<RDM<1>> rdm1_;
    std::shared_ptr<RDM<2>> rdm2_;

    std::shared_ptr<const Coeff> coeff_;

    void print_header() const;
    void print_iteration(int iter, int miter, int tcount, const std::vector<double> energy, const double error, const double time) const;
    void common_init();

    void mute_stdcout();
    void resume_stdcout();

    std::shared_ptr<const Matrix> hcore_;

    std::shared_ptr<const Coeff> update_coeff(const std::shared_ptr<const Matrix> cold, std::shared_ptr<const Matrix> natorb) const;

    // energy
    std::vector<double> energy_;
    double rms_grad_;

    std::shared_ptr<Matrix> Qvec(const int n, const int m, std::shared_ptr<const Matrix> c, const size_t nclosed) const;
    double check_symmetric(std::shared_ptr<Matrix>& mat) const;

  public:
    ASD_OO(std::shared_ptr<const PTree> idat, std::shared_ptr<Dimer> dimer);

    virtual ~ASD_OO();

    virtual void compute() override = 0;

    std::shared_ptr<const Reference> ref() const { return ref_; };
    virtual std::shared_ptr<const Reference> conv_to_ref() const override;

    // functions to retrieve protected members
    int nocc() const { return nocc_; }
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }

    int nactA() const { return nactA_; }
    int nactB() const { return nactB_; }
    std::array<int,3> rasA() const { return rasA_; }
    std::array<int,3> rasB() const { return rasB_; }

    int nvirt() const { return nvirt_; }
    int nbasis() const { return nbasis_; }
    int nstate() const { return nstate_; }
    int max_iter() const { return max_iter_; }
    double thresh() const { return thresh_; }

    double energy(const int i) const { return energy_[i]; };
    double rms_grad() const { return rms_grad_; };

    std::shared_ptr<const Matrix> hcore() const { return hcore_; };
    std::shared_ptr<const Coeff> coeff() const { return coeff_; };
};

}

#endif
