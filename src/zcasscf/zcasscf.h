//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_ZCASSCF_ZCASSCF_H
#define __SRC_ZCASSCF_ZCASSCF_H

#include <src/zfci/zharrison.h>
#include <src/casscf/rotfile.h>
#include <src/wfn/method.h>
#include <src/rel/reloverlap.h>
#include <src/math/bfgs.h>
#include <src/math/step_restrict_bfgs.h>

namespace bagel {

class ZCASSCF : public Method {
  protected:
    int nneg_;
    int nocc_;
    int nclosed_;
    int nact_;
    int nvirt_;
    int nvirtnr_;
    int nbasis_;

    int charge_;

    bool gaunt_;
    bool breit_;

    double thresh_;
    double thresh_micro_;
    std::complex<double> level_shift_;

    int nstate_;

    int max_iter_;
    int max_micro_iter_;

    std::shared_ptr<const ZMatrix> coeff_;
    std::shared_ptr<const Matrix>  nr_coeff_;

    void print_header() const;
    void print_iteration(int iter, int miter, int tcount, const std::vector<double> energy, const double error, const double time) const;

    void init();
    void init_kramers_coeff(std::shared_ptr<const ZMatrix> hcore, std::shared_ptr<const RelOverlap> overlap);

    void mute_stdcout(const bool fci) const;
    void resume_stdcout() const;

    std::shared_ptr<ZHarrison> fci_;
    std::shared_ptr<const ZMatrix> active_fock(std::shared_ptr<const ZMatrix>) const;
    std::shared_ptr<const ZMatrix> transform_rdm1() const;

    // energy
    std::vector<double> energy_;
    double micro_energy_;

    // internal function
    void kramers_adapt(std::shared_ptr<ZRotFile> o, const int nvirt) const;
    void kramers_adapt(std::shared_ptr<ZMatrix> o, const int nvirt) const;

    void zero_positronic_elements(std::shared_ptr<ZRotFile> rot);

    std::shared_ptr<ZMatrix> project_non_rel_coeff(std::shared_ptr<const RelOverlap> overlap, std::shared_ptr<const ZMatrix> hcore, std::shared_ptr<ZMatrix> kramers_coeff) const;
    std::shared_ptr<ZMatrix> nonrel_to_relcoeff(std::shared_ptr<const RelOverlap> overlap, const bool stripes = true) const;

  public:
    ZCASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref = nullptr);

    virtual void compute() override = 0;

    // TODO
    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }

//  private:
    // TODO debug only. All implemented in zcasscf_debug.cc. Will be removed once everything works.
    void ___debug___orbital_rotation(const bool kramers);
    void ___debug___print_gradient(std::shared_ptr<const ZRotFile> grad, const bool kramers) const;
    void ___debug___compute_hessian(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr, const bool kramers) const;
    // returns [x,y] = (xx|yy) (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,y] = (xy|yx) (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,y] = (xx'|y'y) (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb_kramers(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,y] = (xy|y'x') (x is an index of coeffa, and y is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange_kramers(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool off_diagonal = false, const bool diagonal = false) const;

    // returns [x,t,u] = (xx|tu) (x is an index of coeffa, and coeffi should be active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb_active(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [x,t,u] = (xu|tx) (x is an index of coeffa, and coeffi should be active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange_active(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool with_kramers = false) const;
    // returns [s,t,u,v] = (st|uv) (stuv are indicies of coeffi and coeffi should be active)
    std::shared_ptr<ZMatrix> ___debug___all_integrals_coulomb_active(std::shared_ptr<const ZMatrix> coeffi) const;
    // returns [a,t] = (aa|vw) * G(vw,tt) (a is index of coeffa, and t is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_2rdm_contraction_coulomb(std::shared_ptr<const ZMatrix> coeffa) const;
    // returns [a,t] = (aw|va) * G(vw,tt) (a is index of coeffa, and t is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_2rdm_contraction_exchange(std::shared_ptr<const ZMatrix> coeffa, const bool with_kramers = false) const;
    // returns [a,t] = (ka a|vw) * G(vw,kt t) (a is index of coeffa, and t is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_coulomb_active_kramers(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool closed_active = false) const;
    // returns [a,t] = (ka w|v a) * G(vw,kt t) (a is index of coeffa, and t is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_integrals_exchange_active_kramers(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool closed_active = false) const;
    // returns [a,i] = (a a|i u) * D(iu) (a is index of coeffa, and i is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_1rdm_contraction_coulomb(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool with_kramers = false) const;
    // returns [a,i] = (a u|i a) * D(iu) (a is index of coeffa, and i is active)
    std::shared_ptr<ZMatrix> ___debug___diagonal_1rdm_contraction_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi, const bool with_kramers = false) const;
    // returns [a,i] = [ (i a|v ka)  - (i ka|v a) ] * {^{A}D}_{v ki} + [ (ki ka|v a) - (ki a|v ka) ] * {^{A}D}_{v i} (a is index of coeffa, and i is active)
    std::shared_ptr<ZMatrix> ___debug___closed_active_offdiagonal_1rdm_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    std::shared_ptr<ZMatrix> ___debug___closed_active_offdiagonal_2rdm_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    std::shared_ptr<ZMatrix> ___debug___closed_active_diagonal_1rdm_contraction_exchange(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // return FCI energy from transformed 1&2RDM
    double ___debug___recompute_fci_energy(std::shared_ptr<const ZMatrix> cfock) const;

    // CLOSED-ACTIVE debugging
    // returns M(i,t) = G^(1,1)_(ti,ti) (i is an index of coeffi, and t is active)
    std::shared_ptr<ZMatrix> ___debug___closed_active_diagonal_hessian(std::shared_ptr<const ZMatrix> coeffi, std::shared_ptr<const ZMatrix> coefft, std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr, const bool verbose = false) const;
    // returns M(i,t) = G^(1,2)_{kt ki, ti) (i is an index of coeffi, and t is active)
    std::shared_ptr<ZMatrix> ___debug___closed_active_diagonal_hessian_kramers(std::shared_ptr<const ZMatrix> coeffi, std::shared_ptr<const ZMatrix> coefft, const bool verbose) const;
    // returns M(a,i) = G^(1,2)_(kt i, kt i) (i is an index of coeffi, and t is active)
    std::shared_ptr<ZMatrix> ___debug___closed_active_offdiagonal_hessian_kramers(std::shared_ptr<const ZMatrix> coeffi, std::shared_ptr<const ZMatrix> coefft, const bool verbose) const;
    // returns M(a,b) = [ (t u|a b) - (t b|a u) ] D(t u) (ab are indices of coeffa and tu are active)
    std::shared_ptr<ZMatrix> ___debug___active_fock(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> rdm1) const;
    // returns M(t,t) = G(t u,v w) (t u|v w) (tuvw are active)
    std::shared_ptr<ZMatrix> ___debug___active_qvec_byhand(std::shared_ptr<const ZMatrix> coefft) const;
    // returns M(a,i) = (a i|a i) (where a is an index of coeffa and i is an index of coeffi)
    std::shared_ptr<ZMatrix> ___debug___offdiagonal_exchange_integrals(std::shared_ptr<const ZMatrix> coeffa, std::shared_ptr<const ZMatrix> coeffi) const;
    // returns "optimal" level shift
    std::complex<double> find_level_shift(std::shared_ptr<const ZRotFile> rotmat) const;
    // returns energy ratio needed for trust radius updates
    double trust_radius_energy_ratio(const int iter, const std::vector<double> energy, std::shared_ptr<ZRotFile> a, std::shared_ptr<ZRotFile> v, std::shared_ptr<ZRotFile> grad) const;
    // returns matrix of orbital rotation parameters
    std::shared_ptr<ZRotFile> ___debug___microiterations(std::shared_ptr<ZRotFile> xlog, std::shared_ptr<ZRotFile> grad, std::shared_ptr<SRBFGS<ZRotFile>> bfgs, double trust_radius, const int iter) const;
    // function to optimize only the electronic-electronic or electronic-positronic type rotations
    std::shared_ptr<ZRotFile> ___debug___optimize_subspace_rotations(std::vector<double> energy, std::shared_ptr<const ZRotFile> grad, std::shared_ptr<const ZRotFile> rot, std::shared_ptr<SRBFGS<ZRotFile>> srbfgs, bool optimize_electrons = true);
    // function to copy electronic rotations from a rotation file TODO: make lambda
    std::shared_ptr<ZRotFile> ___debug___copy_electronic_rotations(std::shared_ptr<const ZRotFile> rot) const;
    // function to copy positronic rotations from a rotation file TODO: make lambda
    std::shared_ptr<ZRotFile> ___debug___copy_positronic_rotations(std::shared_ptr<const ZRotFile> rot) const;
    // function to compute energy and gradient of a given coeff
    std::shared_ptr<ZRotFile> ___debug___compute_energy_and_gradients(std::shared_ptr<const ZMatrix> coeff, std::shared_ptr<const ZMatrix> hcore);
    // line search function specific to ZCAS
    double  ___debug___line_search(std::shared_ptr<ZRotFile> grad, std::shared_ptr<ZRotFile> delta, std::shared_ptr<const ZMatrix> hcore);
};

}

#endif
