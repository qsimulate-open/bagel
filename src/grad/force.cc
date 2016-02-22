//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: force.cc
// Copyright (C) 2015 Toru Shiozaki
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
#include <src/grad/force.h>
#include <src/grad/gradeval.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

Force::Force(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}


void Force::compute() {
  auto input = idata_->get_child("method");

  shared_ptr<const Reference> ref = ref_;
  auto m = input->begin();
  for ( ; m != --input->end(); ++m) {
    const std::string title = to_lower((*m)->get<std::string>("title", ""));
    if (title != "molecule") {
      shared_ptr<Method> c = construct_method(title, *m, geom_, ref);
      if (!c) throw runtime_error("unknown method in forceimization");
      c->compute();
      ref = c->conv_to_ref();
    } else {
      geom_ = make_shared<const Geometry>(*geom_, *m);
      if (ref) ref = ref->project_coeff(geom_);
    }
  }
  auto cinput = make_shared<PTree>(**m);
  cinput->put("gradient", true);

  const string method = to_lower(cinput->get<string>("title", ""));

  if (method == "uhf") {

    auto force = make_shared<GradEval<UHF>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "rohf") {

    auto force = make_shared<GradEval<ROHF>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "hf") {

    auto force = make_shared<GradEval<RHF>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "ks") {

    auto force = make_shared<GradEval<KS>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "dhf") {

    auto force = make_shared<GradEval<Dirac>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "mp2") {

    auto force = make_shared<GradEval<MP2Grad>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "dmp2") {

    auto force = make_shared<GradEval<DMP2Grad>>(cinput, geom_, ref_);
    force->compute();

  } else if (method == "casscf") {
    string algorithm = cinput->get<string>("algorithm", "");
    // in case of SS-CASSCF
    if (cinput->get<int>("nstate", 1) == 1) {
      if (algorithm == "superci" || algorithm == "") {
        auto force = make_shared<GradEval<SuperCI>>(cinput, geom_, ref_);
        force->compute();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    // in case of SA-CASSCF
    } else {
      if (algorithm == "superci" || algorithm == "") {
        auto force = make_shared<GradEval<SuperCIGrad>>(cinput, geom_, ref_);
        force->compute();
      } else {
        throw runtime_error("unknown CASSCF algorithm specified.");
      }
    }
  } else if (method == "caspt2") {

    auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_);
    force->compute();

  }

}
