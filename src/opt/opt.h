//
// Newint - Parallel electron correlation program.
// Filename: opt.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_OPT_OPT_H
#define __SRC_OPT_OPT_H

#include <map>
#include <memory>
#include <string>
#include <fstream>
#include <chrono>
#include <src/util/bfgs.h>
#include <src/scf/geometry.h>
#include <src/grad/gradfile.h>
#include <src/grad/gradeval.h>

template<typename T>
class Opt {
  protected:
    // entire input
    const std::shared_ptr<const InputData> idata_;
    // options for T
    std::multimap<std::string, std::string> input_;
    std::shared_ptr<Geometry> current_;
    std::shared_ptr<BFGS<GradFile> > bfgs_;

    int iter_;

    // somehow using raw pointers
    std::streambuf* backup_stream_;
    std::ofstream* ofs_;

    // TODO make it adjustable from the input
    const double thresh_;
    static const int maxiter_ = 10;

    // reference geometry
    const std::shared_ptr<const GradFile> refgeom_;

  public:
    Opt(std::shared_ptr<const InputData> idat, std::multimap<std::string, std::string>& inp, const std::shared_ptr<Geometry> geom)
      : idata_(idat), input_(inp), current_(geom), iter_(0), backup_stream_(NULL), thresh_(1.0e-6), refgeom_(new GradFile(geom->xyz())) {
      std::shared_ptr<GradFile> denom(new GradFile(geom->natom(), 1.0));
      bfgs_ = std::shared_ptr<BFGS<GradFile> >(new BFGS<GradFile>(denom));
    };
    ~Opt() {};

    bool next() {
      if (iter_ > 0) mute_stdcout(); 
      auto tp1 = std::chrono::high_resolution_clock::now();
      GradEval<T> eval(input_, current_);
      if (iter_ == 0) {
        print_header();
        mute_stdcout();
      }
      // current geom and grad in the cartesian coordinate
      std::shared_ptr<GradFile> cgrad = eval.compute(); 
      std::shared_ptr<GradFile> cgeom(new GradFile(current_->xyz()));
      std::shared_ptr<GradFile> displ;
      if (true) {
        // x = BX
        // g = (BT)^-1 X
        *cgeom -= *refgeom_;
        std::array<std::unique_ptr<double[]>,2> b;
        b  = current_->compute_internal_coordinate();
        std::shared_ptr<GradFile> dgeom = cgeom->transform(b[0], false);
        std::shared_ptr<GradFile> dgrad = cgrad->transform(b[1], false);
        displ = bfgs_->extrapolate(dgrad, dgeom);

        // self consistent cycle (Eq. 6 of JCP 105,192)
        std::shared_ptr<GradFile> previous = displ->clone();
        for (int i = 0; i != maxiter_; ++i) {
          std::shared_ptr<GradFile> tmp = displ->transform(b[1], true);
          std::shared_ptr<Geometry> tmpgeom = std::shared_ptr<Geometry>(new Geometry(*current_, tmp->xyz(), input_, false, true));
          b  = tmpgeom->compute_internal_coordinate();
          if ((*previous-*tmp).norm() < 1.0e-10) break;
          previous = tmp;
        }
        displ = displ->transform(b[1], true);
      } else { 
        displ = bfgs_->extrapolate(cgrad, cgeom);
      }
      const double gradnorm = cgrad->norm();
      const double disnorm = displ->norm();
      const bool converged = gradnorm < thresh_ && disnorm < thresh_;
      if (!converged) {
        displ->scale(-1.0);
        if (iter_ == 0) displ->scale(0.01);
        current_ = std::shared_ptr<Geometry>(new Geometry(*current_, displ->xyz(), input_));
        current_->print_atoms();
      }

      resume_stdcout(); 
      auto tp2 = std::chrono::high_resolution_clock::now();
      print_iteration(eval.energy(), gradnorm, disnorm, tp1, tp2);
      
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
     
    void print_iteration(const double energy, const double residual, const double step,
                         const std::chrono::high_resolution_clock::time_point tp1,
                         const std::chrono::high_resolution_clock::time_point tp2) const {
      auto dr = std::chrono::duration_cast<std::chrono::milliseconds>(tp2-tp1);
      std::cout << std::setw(8) << iter_ << std::setw(20) << std::setprecision(8) << std::fixed << energy
                                         << std::setw(20) << std::setprecision(8) << std::fixed << residual  
                                         << std::setw(15) << std::setprecision(8) << std::fixed << step  
                                         << std::setw(12) << std::setprecision(2) << std::fixed << dr.count()*0.001 << std::endl;
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

#endif
