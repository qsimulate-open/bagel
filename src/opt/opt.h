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
#include <src/util/bfgs.h>
#include <src/scf/geometry.h>
#include <src/grad/gradient.h>
#include <src/grad/gradeval.h>

template<typename T>
class Opt {
  protected:
    std::multimap<std::string, std::string> input_;
    std::shared_ptr<Geometry> current_;
    std::shared_ptr<BFGS<Gradient> > bfgs_;

  public:
    Opt(std::multimap<std::string, std::string>& idata, const std::shared_ptr<Geometry> geom) : input_(idata), current_(geom) {
      std::shared_ptr<Gradient> denom(new Gradient(geom->natom()*3, 1.0));
      bfgs_ = std::shared_ptr<BFGS<Gradient> >(new BFGS<Gradient>(denom));
    };
    ~Opt() {};

    void next() {
      GradEval<T> opt(input_, current_);
      std::shared_ptr<Gradient> cgrad = opt.compute(); 
    };


};

#endif
