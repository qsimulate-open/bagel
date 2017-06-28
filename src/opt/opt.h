//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: opt.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_OPT_OPT_H
#define __SRC_OPT_OPT_H

#include <functional>
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/wfn/get_energy.h>
#include <src/opt/constraint.h>
#include <src/util/muffle.h>

namespace bagel {

class Opt {
  protected:
    // entire input
    const std::shared_ptr<const PTree> idata_;
    // options for T
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> current_;
    std::shared_ptr<const Reference> prev_ref_;

    int target_state_;
    std::string opttype_;
    int target_state2_;

    int iter_;

    mutable std::shared_ptr<Muffle> muffle_;

    std::string algorithm_;
    std::string method_;
    std::string hess_update_;

    int maxiter_;
    int maxziter_;
    double thresh_grad_;
    double thresh_displ_;
    double thresh_echange_;
    double maxstep_;
    bool scratch_;

    bool numerical_;

    std::array<std::shared_ptr<const Matrix>,3> bmat_;
    std::array<std::shared_ptr<const Matrix>,4> bmat_red_;

    // constraints
    bool constrained_;
    std::vector<std::shared_ptr<const OptConstraint>> constraints_;
    bool explicit_bond_;
    std::vector<std::shared_ptr<const OptExpBonds>> bonds_;

    // whether we use a delocalized internal coordinate or not
    bool internal_;
    bool redundant_;
    int dispsize_;
    // whether we use adaptive stepsize or not
    bool adaptive_;
    // whether we use ab initio hessian or approximate hessian
    bool hess_approx_;
    size_t size_;
    // nonadiabatic coupling type used in conical
    int nacmtype_;
    double thielc3_, thielc4_;
    // MEP direction
    int mep_direction_;

    Timer timer_;

    // some global values needed for quasi-newton optimizations
    double en_;
    double predictedchange_;
    double predictedchange_prev_;

    std::shared_ptr<GradFile> grad_;
    std::shared_ptr<Matrix> hess_;
    std::shared_ptr<XYZFile> displ_;

    // history
    std::vector<double> prev_en_;
    std::vector<std::shared_ptr<const GradFile>> prev_grad_;
    std::vector<std::shared_ptr<const XYZFile>> prev_xyz_;
    std::vector<std::shared_ptr<const XYZFile>> prev_displ_;
    std::shared_ptr<const GradFile> prev_grad_internal_;

    // protected compute module (changes object)
    void compute_optimize();
    void compute_mep(std::shared_ptr<XYZFile> mep_start);

    // const internal functions
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_grad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_grad_energy(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_mecigrad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_mdcigrad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,std::shared_ptr<GradFile>> get_euclidean_dist(std::shared_ptr<const XYZFile> a, std::shared_ptr<const XYZFile> refgeom) const;

    std::tuple<double,double,std::shared_ptr<XYZFile>> get_step() const;
    std::shared_ptr<XYZFile> get_step_nr() const;
    std::shared_ptr<XYZFile> get_step_ef() const;
    std::shared_ptr<XYZFile> get_step_ef_pn() const;
    std::tuple<double,double,std::shared_ptr<XYZFile>> get_step_rfo() const;

    std::shared_ptr<XYZFile> iterate_displ() const;

    std::shared_ptr<Matrix> hessian_update() const;
    std::shared_ptr<Matrix> hessian_update_bfgs(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> hs) const;
    std::shared_ptr<Matrix> hessian_update_sr1(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> z) const;
    std::shared_ptr<Matrix> hessian_update_psb(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> z) const;

    double do_adaptive() const;

  public:
    Opt(std::shared_ptr<const PTree> idat, std::shared_ptr<const PTree> inp, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref);

    ~Opt() {
      print_footer();
      current_->print_atoms();
    }

    void compute();

    void print_header() const;
    void print_footer() const { std::cout << std::endl << std::endl; }

    void print_iteration_energy(const double residual, const double time) const;
    void print_iteration_conical(const double residual, const double param, const double time) const;

    void print_history_molden() const;

    void print_iteration(const double residual, const double param, const double time) const;

    std::shared_ptr<const Geometry> geometry() const { return current_; }
    std::shared_ptr<const Reference> conv_to_ref() const { return prev_ref_; }

};

}
#endif
