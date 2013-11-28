//
// BAGEL - Parallel electron correlation program.
// Filename: optimize.cc
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

#include <string>
#include <src/opt/optimize.h>
#include <src/opt/opt.h>

using namespace std;
using namespace bagel;

Optimize::Optimize(const shared_ptr<const PTree> idata, shared_ptr<const Geometry> g) : idata_(idata), geom_(g) {
  maxiter_ = idata->get<int>("maxiter", 100);

}


void Optimize::compute() {
  auto lastmethod = *idata_->get_child("method")->rbegin();
  string method = lastmethod->get<string>("title", "");
  transform(method.begin(), method.end(), method.begin(), ::tolower);
  if (method.empty())
    throw runtime_error("title is missing in one of the input blocks (opt)");

  shared_ptr<const PTree> methodblock = idata_->get_child("method");

  if (method == "uhf") {

    auto opt = make_shared<Opt<UHF>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "rohf") {

    auto opt = make_shared<Opt<ROHF>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "hf") {

    auto opt = make_shared<Opt<SCF>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "ks") {

    auto opt = make_shared<Opt<KS>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "dhf") {

    auto opt = make_shared<Opt<Dirac>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "mp2") {

    auto opt = make_shared<Opt<MP2Grad>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "dmp2") {

    auto opt = make_shared<Opt<DMP2Grad>>(idata_, methodblock, geom_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "casscf") {
    string algorithm = lastmethod->get<string>("algorithm", "");
    // in case of SS-CASSCF
    if (lastmethod->get<int>("nstate", 1) == 1) {
      if (algorithm == "superci" || algorithm == "") {
        auto opt = make_shared<Opt<SuperCI>>(idata_, methodblock, geom_);
        opt->compute();
        geom_ = opt->geometry();
      } else if (algorithm == "werner" || algorithm == "knowles") {
        auto opt = make_shared<Opt<WernerKnowles>>(idata_, methodblock, geom_);
        opt->compute();
        geom_ = opt->geometry();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    // in case of SA-CASSCF
    } else {
      if (algorithm == "superci" || algorithm == "") {
        auto opt = make_shared<Opt<SuperCIGrad>>(idata_, methodblock, geom_);
        opt->compute();
        geom_ = opt->geometry();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    }
  }

}
