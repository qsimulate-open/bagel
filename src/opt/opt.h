//
// BAGEL - Parallel electron correlation program.
// Filename: opt.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_OPT_OPT_H
#define __SRC_OPT_OPT_H

#include <fstream>
#include <src/util/bfgs.h>
#include <src/util/timer.h>
#include <src/grad/gradeval.h>

namespace bagel {

template<typename T>
class Opt {
  protected:
    // entire input
    const std::shared_ptr<const PTree> idata_;
    // options for T
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> current_;
    std::shared_ptr<BFGS<GradFile>> bfgs_;

    int iter_;

    // somehow using raw pointers
    std::streambuf* backup_stream_;
    std::ofstream* ofs_;

    // TODO make it adjustable from the input
    const double thresh_;

    static const int maxiter_ = 10;
    static const bool nodf = true;
    static const bool rotate = false;

    // reference geometry
    const std::shared_ptr<const GradFile> refgeom_;

    std::array<std::unique_ptr<double[]>,2> bmat_;

    // whether we use a delocalized internal coordinate or not
    bool internal_;

  public:
    Opt(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const PTree> inp, const std::shared_ptr<const Geometry> geom)
      : idata_(idat), input_(inp), current_(geom), iter_(0), backup_stream_(nullptr), thresh_(1.0e-5), refgeom_(std::make_shared<GradFile>(geom->xyz())) {
      bfgs_ = std::make_shared<BFGS<GradFile>>(std::make_shared<const GradFile>(geom->natom(), 1.0));
      bmat_ = current_->compute_internal_coordinate();

      internal_ = inp->get<bool>("internal", true);
    }

    bool next() {
      if (iter_ > 0) mute_stdcout();
      Timer timer;
      GradEval<T> eval(input_, current_);
      if (iter_ == 0) {
        print_header();
        mute_stdcout();
      }
      // current geom and grad in the cartesian coordinate
      std::shared_ptr<const GradFile> cgrad = eval.compute();
      std::shared_ptr<const GradFile> cgeom = std::make_shared<GradFile>(current_->xyz());
      std::shared_ptr<GradFile> displ;
      if (internal_) {
        std::shared_ptr<const GradFile> dgeom = cgeom->transform(bmat_[0], false);
        std::shared_ptr<const GradFile> dgrad = cgrad->transform(bmat_[1], false);
        displ = bfgs_->extrapolate(dgrad, dgeom);
        // TODO  I haven't understood why this update (internal displacement to cartesian) can be iterative!!
#if 0
        // self consistent cycle (Eq. 6 of JCP 105,192)
        std::shared_ptr<GradFile> previous = displ;
        for (int i = 0; i != maxiter_; ++i) {
        }
#else
        displ = displ->transform(bmat_[1], true);
#endif
      } else {
        displ = bfgs_->extrapolate(cgrad, cgeom);
      }
      const double gradnorm = cgrad->norm();
      const double disnorm = displ->norm();
      const bool converged = gradnorm < thresh_ && disnorm < thresh_;
      if (!converged) {
        displ->scale(-1.0);
        if (iter_ == 0) displ->scale(0.01);
        current_ = std::make_shared<const Geometry>(*current_, displ->xyz(), input_);
        current_->print_atoms();
      }

      resume_stdcout();
      print_iteration(eval.energy(), gradnorm, disnorm, timer.tick());

      ++iter_;
      if (converged) { print_footer(); current_->print_atoms(); }
      return gradnorm < thresh_ && disnorm < thresh_;
    };

    void print_footer() const { std::cout << std::endl << std::endl; };
    void print_header() const {
        std::cout << std::endl << "  *** Geometry optimization started ***" << std::endl <<
                                  "     iter           energy             res norm      step norm"
        << std::endl << std::endl;
    };

    void print_iteration(const double energy, const double residual, const double step, const double time) const {
      std::cout << std::setw(8) << iter_ << std::setw(20) << std::setprecision(8) << std::fixed << energy
                                         << std::setw(20) << std::setprecision(8) << std::fixed << residual
                                         << std::setw(15) << std::setprecision(8) << std::fixed << step
                                         << std::setw(12) << std::setprecision(2) << std::fixed << time << std::endl;
    };

    void mute_stdcout() {
      ofs_ = new std::ofstream("opt.log",(backup_stream_ ? std::ios::app : std::ios::trunc));
      backup_stream_ = std::cout.rdbuf(ofs_->rdbuf());
    };

    void resume_stdcout() {
      std::cout.rdbuf(backup_stream_);
      delete ofs_;
    };

    std::shared_ptr<const Geometry> geometry() const { return current_; };

};

}

#endif
