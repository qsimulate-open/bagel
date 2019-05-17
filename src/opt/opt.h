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
#include <src/opt/qmmm.h>
#include <src/opt/optinfo.h>

namespace bagel {

class Opt {
  protected:
    // input and output
    // entire input
    const std::shared_ptr<const PTree> idata_;

    // options for method block
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> current_;
    std::shared_ptr<const Reference> prev_ref_;
    std::string method_;

    // input parameters kept constant
    std::shared_ptr<const OptInfo> optinfo_;

    // QM/MM
    std::shared_ptr<const QMMM> qmmm_driver_;

    // output
    mutable std::shared_ptr<Muffle> muffle_;
    Timer timer_;

    // global values that change during the optimization
    // internal - Cartesian transformation
    std::array<std::shared_ptr<const Matrix>,3> bmat_;
    std::vector<std::vector<int>> bondlist_;
    std::array<std::shared_ptr<const Matrix>,5> bmat_red_;

    // size of internals
    int iter_;
    int dispsize_;
    size_t size_;

    // some values for Newton-Raphson optimization
    double maxstep_;
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
    std::vector<std::shared_ptr<const XYZFile>> prev_xyz_internal_;
    std::vector<std::shared_ptr<const XYZFile>> prev_displ_;
    std::vector<std::shared_ptr<const GradFile>> prev_grad_internal_;

    // protected compute module (changes object)
    void compute_optimize();
    void compute_mep(std::shared_ptr<XYZFile> mep_start);

    // const internal functions
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_grad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_grad_energy(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_mecigrad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,double,std::shared_ptr<const Reference>,std::shared_ptr<GradFile>> get_mdcigrad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref) const;
    std::tuple<double,std::shared_ptr<GradFile>> get_euclidean_dist(std::shared_ptr<const XYZFile> a, std::shared_ptr<const XYZFile> refgeom) const;
    std::tuple<std::shared_ptr<PTree>,std::shared_ptr<const Reference>,std::shared_ptr<const Geometry>> get_grad_input() const;

    std::tuple<double,double,std::shared_ptr<XYZFile>> get_step() const;
    std::shared_ptr<XYZFile> get_step_nr() const;
    std::shared_ptr<XYZFile> get_step_ef() const;
    std::shared_ptr<XYZFile> get_step_ef_pn(const int mode) const;
    std::tuple<double,double,std::shared_ptr<XYZFile>> get_step_rfo() const;

    std::shared_ptr<XYZFile> iterate_displ() const;

    std::shared_ptr<Matrix> hessian_update() const;
    std::shared_ptr<Matrix> hessian_update_bfgs(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> hs) const;
    std::shared_ptr<Matrix> hessian_update_sr1(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> z) const;
    std::shared_ptr<Matrix> hessian_update_psb(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> z) const;
    std::shared_ptr<Matrix> hessian_update_bofill(std::shared_ptr<const GradFile> y, std::shared_ptr<const GradFile> s, std::shared_ptr<const GradFile> z) const;

    double do_adaptive(const int iter) const;

  public:
    Opt(std::shared_ptr<const PTree> idat, std::shared_ptr<const PTree> inp, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref);

    ~Opt() {
      print_footer();
      current_->print_atoms();
    }

    void compute();

    void print_header() const;
    void print_footer() const { std::cout << std::endl << std::endl; }

    void print_iteration_energy(const int iter, const double residual, const double time) const;
    void print_iteration_conical(const int iter, const double residual, const double param, const double time) const;

    void print_history_molden() const;

    void print_iteration(const int iter, const double residual, const double param, const double time) const;

    std::shared_ptr<const OptInfo> optinfo() const { return optinfo_; }
    std::shared_ptr<const Geometry> geometry() const { return current_; }
    std::shared_ptr<const Reference> conv_to_ref() const { return prev_ref_; }

};

}
#endif
