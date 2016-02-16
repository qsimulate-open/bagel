//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: optimize.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <string>
#include <src/opt/optimize.h>
#include <src/opt/opt.h>

using namespace std;
using namespace bagel;

Optimize::Optimize(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {
  maxiter_ = idata->get<int>("maxiter", 100);

}


void Optimize::compute() {
  auto lastmethod = *idata_->get_child("method")->rbegin();
  const string method = to_lower(lastmethod->get<string>("title", ""));
  if (method.empty())
    throw runtime_error("title is missing in one of the input blocks (opt)");

  shared_ptr<const PTree> methodblock = idata_->get_child("method");

  if (method == "uhf") {

    auto opt = make_shared<Opt<UHF>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "rohf") {

    auto opt = make_shared<Opt<ROHF>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "hf") {

    auto opt = make_shared<Opt<RHF>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "ks") {

    auto opt = make_shared<Opt<KS>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "dhf") {

    auto opt = make_shared<Opt<Dirac>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "mp2") {

    auto opt = make_shared<Opt<MP2Grad>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "dmp2") {

    auto opt = make_shared<Opt<DMP2Grad>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  } else if (method == "casscf") {
    string algorithm = lastmethod->get<string>("algorithm", "");
    // in case of SS-CASSCF
    if (lastmethod->get<int>("nstate", 1) == 1) {
      if (algorithm == "superci" || algorithm == "") {
        auto opt = make_shared<Opt<SuperCI>>(idata_, methodblock, geom_, ref_);
        opt->compute();
        geom_ = opt->geometry();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    // in case of SA-CASSCF
    } else {
      if (algorithm == "superci" || algorithm == "") {
        auto opt = make_shared<Opt<SuperCIGrad>>(idata_, methodblock, geom_, ref_);
        opt->compute();
        geom_ = opt->geometry();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    }
  } else if (method == "caspt2") {

    auto opt = make_shared<Opt<CASPT2Grad>>(idata_, methodblock, geom_, ref_);
    opt->compute();
    geom_ = opt->geometry();

  }

}
