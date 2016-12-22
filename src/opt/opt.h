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
#include <src/wfn/construct_method.h>
#include <src/alglib/optimization.h>

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

    // somehow using raw pointers
    std::streambuf* backup_stream_;
    std::ofstream* ofs_;

    std::string algorithm_;
    std::string method_;

    int maxiter_;
    double thresh_;
    double maxstep_;
    bool scratch_;

    std::array<std::shared_ptr<const Matrix>,3> bmat_;
    std::array<std::shared_ptr<const Matrix>,4> bmat_red_;

    // whether we use alglib or not
    bool alglib_;
    // whether we use a delocalized internal coordinate or not
    bool internal_;
    bool redundant_;
    int dispsize_;
    // whether we use adaptive stepsize or not
    bool adaptive_;
    // whether we save reference or not
    bool refsave_;
    std::string refname_;
    size_t size_;
    // nonadiabatic coupling type used in conical
    int nacmtype_;
    double thielc3_, thielc4_;

    Timer timer_;

    // some global values needed for ALGLIB-based optimizations
    void evaluate_alglib(const alglib::real_1d_array& x, double& en, alglib::real_1d_array& grad, void* ptr);
    using eval_type = std::function<void(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*)>;

    // some global values needed for quasi-newton optimizations
    double en_;
    double en_prev_;
    double egap_;
    double predictedchange_;
    double predictedchange_prev_;
    double realchange_;

    std::shared_ptr<GradFile> grad_;
    std::shared_ptr<GradFile> prev_grad_;
    std::shared_ptr<Matrix> hess_;
    std::shared_ptr<XYZFile> displ_;

    // some internal functions 
    std::shared_ptr<GradFile> get_grad(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref);
    std::shared_ptr<GradFile> get_grad_energy(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref);
    std::shared_ptr<GradFile> get_cigrad_bearpark(std::shared_ptr<PTree> cinput, std::shared_ptr<const Reference> ref);

    std::shared_ptr<XYZFile> get_step();
    std::shared_ptr<XYZFile> get_step_steep();
    std::shared_ptr<XYZFile> get_step_nr();
    std::shared_ptr<XYZFile> get_step_rfo();
    std::shared_ptr<XYZFile> get_step_rfos();

    void hessian_update();
    void hessian_update_bfgs(std::shared_ptr<GradFile> y, std::shared_ptr<GradFile> s, std::shared_ptr<GradFile> hs);
    void hessian_update_sr1(std::shared_ptr<GradFile> y, std::shared_ptr<GradFile> s, std::shared_ptr<GradFile> z);
    void hessian_update_psb(std::shared_ptr<GradFile> y, std::shared_ptr<GradFile> s, std::shared_ptr<GradFile> z);

    void do_adaptive();

  public:
    Opt(std::shared_ptr<const PTree> idat, std::shared_ptr<const PTree> inp, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref)
      : idata_(idat), input_(inp), current_(geom), prev_ref_(ref), iter_(0), backup_stream_(nullptr) {

      auto lastmethod = *idat->get_child("method")->rbegin();
      method_ = to_lower(lastmethod->get<std::string>("title", ""));

      target_state_ = idat->get<int>("target", 0);
      internal_ = idat->get<bool>("internal", true);
      redundant_ = idat->get<bool>("redundant", false);
      maxiter_ = idat->get<int>("maxiter", 100);
      maxstep_ = idat->get<double>("maxstep", 0.1);
      scratch_ = idat->get<bool>("scratch", false);
      if (internal_) {
        if (redundant_) 
          bmat_red_ = current_->compute_redundant_coordinate();
        else
          bmat_ = current_->compute_internal_coordinate();
      }
      thresh_ = idat->get<double>("thresh", 5.0e-5);
      algorithm_ = idat->get<std::string>("algorithm", "rfo");
      adaptive_ = idat->get<bool>("adaptive", true);
      opttype_ = idat->get<std::string>("opttype", "energy");
#ifndef DISABLE_SERIALIZATION
      refsave_ = idat->get<bool>("save_ref", false);
      if (refsave_) refname_ = idat->get<std::string>("ref_out", "reference");
#endif

      // For LBFGS or CG optimizer, ALGLIB is true. Otherwise, ALGLIB is false
      if (algorithm_ == "lbfgs" || algorithm_ == "cg")
        alglib_ = true;
      else 
        alglib_ = false;

      if (opttype_ == "conical") {
        target_state2_ = idat->get<int>("target2", 1);
        if (target_state2_ > target_state_) {
          int tmpstate = target_state_;
          target_state_ = target_state2_;
          target_state2_ = tmpstate;
        }
        nacmtype_ = idat->get<int>("nacmtype", 0);
        thielc3_  = idat->get<double>("thielc3", 2.0);
        thielc4_  = idat->get<double>("thielc4", 0.5);
        adaptive_ = false;        // we cannot use it for conical intersection optimization because we do not have a target function
      }
      else if (opttype_ == "transition")
        throw std::runtime_error("We cannot do saddle point optimization now, wait for Hessian coming up...");
      else if (opttype_ != "energy")
        throw std::runtime_error("Optimization type should be: \"energy\", \"transition\" or \"conical\"");
    }

    ~Opt() {
      print_footer();
      current_->print_atoms();
    }

    void compute_alglib();
    void compute_noalglib();
    void compute() {
      if (alglib_)
        compute_alglib();
      else
        compute_noalglib();
    }

    void print_footer() const { std::cout << std::endl << std::endl; };
    void print_header() const {
      if (opttype_ == "energy") {
        std::cout << std::endl << "  *** Geometry optimization started ***" << std::endl <<
                                  "     iter         energy               grad rms       time"
        << std::endl << std::endl;
      } else if (opttype_ == "conical") {
        std::cout << std::endl << "  *** Conical intersection optimization started ***" << std::endl <<
                                  "     iter         energy             gap energy            grad rms       time"
        << std::endl << std::endl;
      }
    }

    void print_iteration_energy(const double residual, const double time) const {
      std::cout << std::setw(7) << iter_ << std::setw(20) << std::setprecision(8) << std::fixed << en_
                                         << std::setw(20) << std::setprecision(8) << std::fixed << residual
                                         << std::setw(12) << std::setprecision(2) << std::fixed << time << std::endl;
    }

    void print_iteration_conical(const double residual, const double time) const {
      std::cout << std::setw(7) << iter_ << std::setw(20) << std::setprecision(8) << std::fixed << en_
                                         << std::setw(20) << std::setprecision(8) << std::fixed << egap_
                                         << std::setw(20) << std::setprecision(8) << std::fixed << residual
                                         << std::setw(12) << std::setprecision(2) << std::fixed << time << std::endl;
    }

    void print_iteration(const double residual, const double time) {
      if (opttype_ == "energy") print_iteration_energy (residual, time);
      else if (opttype_ == "conical") print_iteration_conical (residual, time);
    }

    void mute_stdcout() {
      if (mpi__->rank() == 0) {
        ofs_ = new std::ofstream("opt.log",(backup_stream_ ? std::ios::app : std::ios::trunc));
        backup_stream_ = std::cout.rdbuf(ofs_->rdbuf());
      }
    }

    void resume_stdcout() {
      if (mpi__->rank() == 0) {
        std::cout.rdbuf(backup_stream_);
        delete ofs_;
      }
    }
    std::shared_ptr<const Geometry> geometry() const { return current_; }

};

}
#endif
