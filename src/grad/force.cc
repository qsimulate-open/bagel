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
#include <src/grad/finite.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

Force::Force(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}


shared_ptr<GradFile> Force::compute() {
  auto input = idata_->get_child("method");
  const int target = idata_->get<int>("target", 0);
  const string jobtitle = to_lower(idata_->get<string>("title", ""));
  const int target2= idata_->get<int>("target2", 1);
  const int nacmtype = idata_->get<int>("nacmtype", 0);
  const bool export_grad = idata_->get<bool>("export", false);
  const bool export_single = idata_->get<bool>("export_single", false);
#ifndef DISABLE_SERIALIZATION
  const bool refsave = idata_->get<bool>("save_ref", false);
  string refname;
  if (refsave) refname = idata_->get<string>("ref_out", "reference");
#endif

  shared_ptr<const Reference> ref = ref_;
  auto m = input->begin();
  for ( ; m != --input->end(); ++m) {
    const std::string title = to_lower((*m)->get<std::string>("title", ""));
    if (title != "molecule") {
      shared_ptr<Method> c = construct_method(title, *m, geom_, ref);
      if (!c) throw runtime_error("unknown method in force");
      c->compute();
      ref = c->conv_to_ref();
    } else {
      geom_ = make_shared<const Geometry>(*geom_, *m);
      if (ref) ref = ref->project_coeff(geom_);
    }
  }
  auto cinput = make_shared<PTree>(**m);
  cinput->put("gradient", true);

  numerical_ = idata_->get<bool>("numerical", false);
  if (geom_->dkh()) {
    cout << "  Analytcal gradient not supplied for DKH Hamiltonian." << endl;
    numerical_ = true;
  }
  if (numerical_)
    cout << "  The gradients will be computed with finite difference." << endl;
  else
    cout << "  The gradients will be computed analytically." << endl;

  shared_ptr<GradFile> out;

  const string method = to_lower(cinput->get<string>("title", ""));

  if (!numerical_) {
    if (method == "uhf") {

      auto force = make_shared<GradEval<UHF>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "rohf") {

      auto force = make_shared<GradEval<ROHF>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "hf") {

      auto force = make_shared<GradEval<RHF>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "ks") {

      auto force = make_shared<GradEval<KS>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "dhf") {

      auto force = make_shared<GradEval<Dirac>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();

    } else if (method == "mp2") {

      auto force = make_shared<GradEval<MP2Grad>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "casscf" && jobtitle == "nacme") {

      auto force = make_shared<NacmEval<CASSCF>>(cinput, geom_, ref_, target, target2, nacmtype);
      out = force->compute();
      ref = force->ref();

    } else if (method == "casscf" && jobtitle == "dgrad") {

      auto force = make_shared<DgradEval<CASSCF>>(cinput, geom_, ref_, target, target2);
      out = force->compute();
      ref = force->ref();

    } else if (method == "casscf") {

      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "caspt2" && jobtitle == "nacme") {

      auto force = make_shared<NacmEval<CASPT2Nacm>>(cinput, geom_, ref_, target, target2, nacmtype);
      out = force->compute();
      ref = force->ref();

    } else if (method == "caspt2") {

      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_, target);
      out = force->compute();
      ref = force->ref();
      force_dipole_ = force->dipole();

    }
    else {

        numerical_ = true;
        cout << "  It seems like no analytical gradient method available; moving to finite difference " << endl;

    }
  }

  if (numerical_) {

    const double dx = idata_->get<double>("diffsize", 1.0e-3);

    if (jobtitle == "force") {

      auto force = make_shared<FiniteGrad>(method, cinput, geom_, ref_, target, dx);
      out = force->compute();
      ref = force->ref();

    } else if (jobtitle == "nacme") {

      if (method == "casscf") {

        auto force = make_shared<FiniteNacm<CASSCF>>(method, cinput, geom_, ref_, target, target2, dx);
        out = force->compute();
        ref = force->ref();

      } else if (method == "caspt2") {

        auto force = make_shared<FiniteNacm<CASPT2Energy>>(method, cinput, geom_, ref_, target, target2, dx);
        out = force->compute();
        ref = force->ref();

      }
    }
  }
#ifndef DISABLE_SERIALIZATION
  if (refsave) {
    OArchive archive(refname);
    archive << ref;
  }
#endif

  if (export_grad) {
    if (export_single) {
      shared_ptr<Muffle> wholemuffle;
      wholemuffle = make_shared<Muffle>("FORCE.out");
      wholemuffle->mute();
      cout << setw(20) << setprecision(10) << ref->energy(target) << endl;
      out->print_export();
      wholemuffle->unmute();
    }
    shared_ptr<Muffle> gradmuffle;
    string mufflename = to_upper(jobtitle) + "_" + to_string(target);
    if (jobtitle == "nacme") mufflename += ("_" + to_string(target2));
    mufflename += ".out";
    gradmuffle = make_shared<Muffle>(mufflename);
    gradmuffle->mute();
    out->print_export();
    gradmuffle->unmute();
    shared_ptr<Muffle> enermuffle;
    enermuffle = make_shared<Muffle>("ENERGY.out");
    enermuffle->mute();
    for (int istate = 0; istate != ref->nstate(); ++istate)
      cout << setw(20) << setprecision(10) << ref->energy(istate) << endl;
    enermuffle->unmute();
  }
  return out;
}
