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


#ifndef __SRC_OPT_OPT_H
#define __SRC_OPT_OPT_H

#include <functional> 
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/util/timer.h>
#include <src/grad/gradeval.h>
#include <src/io/moldenout.h>
#include <src/wfn/construct_method.h>
#include <src/alglib/optimization.h>

namespace bagel {

template<typename T>
class Opt {
  protected:
    // entire input
    const std::shared_ptr<const PTree> idata_;
    // options for T
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Geometry> current_;
    std::shared_ptr<const Reference> prev_ref_;

    int iter_;

    // somehow using raw pointers
    std::streambuf* backup_stream_;
    std::ofstream* ofs_;

    std::string algorithm_;

    int maxiter_;
    double thresh_;
    double maxstep_;
    static const bool nodf = true;
    static const bool rotate = false;

    std::array<std::shared_ptr<const Matrix>,2> bmat_;

    // whether we use a delocalized internal coordinate or not
    bool internal_;
    size_t size_;

    Timer timer_;

    void evaluate(const alglib::real_1d_array& x, double& en, alglib::real_1d_array& grad, void* ptr);
    using eval_type = std::function<void(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*)>;

  public:
    Opt(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const PTree> inp, const std::shared_ptr<const Geometry> geom)
      : idata_(idat), input_(inp), current_(geom), iter_(0), backup_stream_(nullptr) {

      internal_ = idat->get<bool>("internal", true);
      maxiter_ = idat->get<int>("maxiter", 100);
      maxstep_ = idat->get<double>("maxstep", 0.1);
      if (internal_)
        bmat_ = current_->compute_internal_coordinate();
      thresh_ = idat->get<double>("thresh", 1.0e-5);
      algorithm_ = idat->get<std::string>("algorithm", "lbfgs");
    }

    ~Opt() {
      print_footer();
      current_->print_atoms();
    }

    void compute() {
      assert(typeid(double) == typeid(double));
      std::shared_ptr<const XYZFile> displ = current_->xyz();
      size_ = internal_ ? bmat_[0]->mdim() : displ->size();

      if (internal_)
        displ = displ->transform(bmat_[0], true);

      try {
        alglib::real_1d_array x;
        x.setcontent(size_, displ->data()); 
        eval_type eval = std::bind(&Opt<T>::evaluate, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

        if (algorithm_ == "cg") {
          alglib::mincgstate state;
          alglib::mincgreport rep;

          alglib::mincgcreate(x, state); 
          alglib::mincgsetcond(state, thresh_*std::sqrt(size_), 0.0, 0.0, maxiter_);
          alglib::mincgsetstpmax(state, maxstep_);

          alglib::mincgoptimize(state, eval);
        } else if (algorithm_ == "lbfgs") {
          alglib::minlbfgsstate state;
          alglib::minlbfgsreport rep;

          alglib::minlbfgscreate(1, x, state); 
          alglib::minlbfgssetcond(state, thresh_*std::sqrt(size_), 0.0, 0.0, maxiter_);
          alglib::minlbfgssetstpmax(state, maxstep_);

          alglib::minlbfgsoptimize(state, eval);
        } else {
          throw std::runtime_error("geometry optimization is implemented only with \"cg\" and \"lbfgs\"");
        }
      } catch (alglib::ap_error e) {
        std::cout << e.msg << std::endl;
        throw std::runtime_error("optimization failed");
      }
    }

    void print_footer() const { std::cout << std::endl << std::endl; };
    void print_header() const {
        std::cout << std::endl << "  *** Geometry optimization started ***" << std::endl <<
                                  "     iter         energy               grad rms       time"
        << std::endl << std::endl;
    }

    void print_iteration(const double energy, const double residual, const double time) const {
      std::cout << std::setw(7) << iter_ << std::setw(20) << std::setprecision(8) << std::fixed << energy
                                         << std::setw(20) << std::setprecision(8) << std::fixed << residual
                                         << std::setw(12) << std::setprecision(2) << std::fixed << time << std::endl;
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


template<class T>
void Opt<T>::evaluate(const alglib::real_1d_array& x, double& en, alglib::real_1d_array& grad, void* ptr) {
  std::shared_ptr<const XYZFile> xyz = current_->xyz();

  // first convert x to the geometry
  auto displ = std::make_shared<XYZFile>(current_->natom());
  assert(size_ == x.length());
  for (int i = 0; i != size_; ++i)
    displ->element(i%3, i/3) = x[i]; // TODO this is just a hack

  if (internal_)
    displ = displ->transform(bmat_[1], false);

  *displ -= *xyz;

  if (iter_ > 0) mute_stdcout();

  // current Geometry
  if (iter_ > 0) {
    current_ = std::make_shared<Geometry>(*current_, displ, std::make_shared<const PTree>()); 
    current_->print_atoms();
    if (internal_)
      bmat_ = current_->compute_internal_coordinate(bmat_[0]);
  }

  // first calculate reference (if needed)
  std::shared_ptr<PTree> cinput; 
  std::shared_ptr<const Reference> ref;
  if (iter_ == 0) {
    auto m = input_->begin();
    for ( ; m != --input_->end(); ++m) {
      std::string title = (*m)->get<std::string>("title", ""); 
      std::transform(title.begin(), title.end(), title.begin(), ::tolower);
      if (title != "molecule") {
        std::shared_ptr<Method> c = construct_method(title, *m, current_, ref);
        if (!c) throw std::runtime_error("unknown method in optimization");
        c->compute();
        ref = c->conv_to_ref();
      } else {
        current_ = std::make_shared<const Geometry>(*current_, *m); 
        if (ref) ref = ref->project_coeff(current_);
      }
    }
    cinput = std::make_shared<PTree>(**m);
  } else {
    ref = prev_ref_->project_coeff(current_);
    cinput = std::make_shared<PTree>(**input_->rbegin());
  }
  cinput->put("gradient", true);

  // then calculate gradients
  GradEval<T> eval(cinput, current_, ref);
  if (iter_ == 0) {
    print_header();
    mute_stdcout();
  }
  // current geom and grad in the cartesian coordinate
  std::shared_ptr<const GradFile> cgrad = eval.compute();
  if (internal_)
    cgrad = cgrad->transform(bmat_[1], true);

  assert(size_ == grad.length());
  for (int i = 0; i != size_; ++i)
    grad[i] = cgrad->element(i%3, i/3); // TODO this is just a hack

  resume_stdcout();

  prev_ref_ = eval.ref();
  ++iter_;
  // returns energy
  en = eval.energy(); 

  // current geometry in a molden file
  MoldenOut mfs("opt.molden");
  mfs << current_;
  print_iteration(en, cgrad->rms(), timer_.tick());
}


}

#endif
