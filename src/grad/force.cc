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
#include <src/wfn/get_energy.h>

using namespace std;
using namespace bagel;

Force::Force(shared_ptr<const PTree> idata, shared_ptr<const Geometry> g, shared_ptr<const Reference> r) : idata_(idata), geom_(g), ref_(r) {

}


shared_ptr<GradFile> Force::compute() {
  string firsttitle = to_lower(idata_->get<string>("title", ""));
  const string jobtitle = (firsttitle == "hessian") ? "force" : firsttitle;
  auto input = idata_->get_child("method");

  shared_ptr<const Reference> ref = ref_;
  auto m = input->begin();
  for ( ; m != --input->end(); ++m) {
    const string title = to_lower((*m)->get<string>("title", ""));
    if (title != "molecule") {
      tie(ignore, ref) = get_energy(title, *m, geom_, ref);
    } else {
      geom_ = make_shared<const Geometry>(*geom_, *m);
      if (ref) ref = ref->project_coeff(geom_);
    }
  }
  auto cinput = make_shared<PTree>(**m);
  cinput->put("_gradient", true);

  numerical_ = idata_->get<bool>("numerical", false);
  if (numerical_)
    cout << "  The gradients will be computed with finite difference." << endl;
  else
    cout << "  The gradients will be computed analytically." << endl;

  shared_ptr<GradFile> out;

  const string method = to_lower(cinput->get<string>("title", ""));
  vector<double> energyvec;

  const bool export_grad = idata_->get<bool>("export", false);
  const bool export_single = idata_->get<bool>("export_single", false);
  const bool compute_dipole = idata_->get<bool>("dipole", false);

  if (jobtitle == "forces") {

    // list of gradients (and NACMEs) to be evaluated. maxziter, nacmtype, target, target2 can be specified
    auto joblist = idata_->get_child("grads");

    if (method == "casscf") {

      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      if (compute_dipole)
        force->compute_dipole();
      for (auto& m : *joblist) {
        const string mtitle= to_lower(m->get<string>("title", "force"));
        auto gradinfo = make_shared<const GradInfo>(m);
        out = force->compute(mtitle, gradinfo);

        if (export_grad)
          force_export(mtitle, gradinfo, energyvec, out, export_single);
      }

    } else if (method == "caspt2") {

      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      if (compute_dipole)
        force->compute_dipole();
      for (auto& m : *joblist) {
        const string mtitle= to_lower(m->get<string>("title", "force"));
        auto gradinfo = make_shared<const GradInfo>(m);
        out = force->compute(mtitle, gradinfo);

        if (export_grad)
          force_export(mtitle, gradinfo, energyvec, out, export_single);
      }

    } else {

      throw runtime_error("multiple force calculations can be only done with CASSCF and CASPT2");

    }

  } else if (!numerical_) {

    auto gradinfo = make_shared<const GradInfo>(idata_);

    if (method == "uhf") {

      auto force = make_shared<GradEval<UHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else if (method == "rohf") {

      auto force = make_shared<GradEval<ROHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else if (method == "hf") {

      auto force = make_shared<GradEval<RHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else if (method == "ks") {

      auto force = make_shared<GradEval<KS>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else if (method == "dhf") {

      auto force = make_shared<GradEval<Dirac>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      out = force->compute(jobtitle, gradinfo);

    } else if (method == "mp2") {

      auto force = make_shared<GradEval<MP2Grad>>(cinput, geom_, ref_);
      out = force->compute(jobtitle, gradinfo);
      ref = force->ref();
      energyvec = force->energyvec();
      force_dipole_ = force->dipole();

    } else if (method == "casscf") {

      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      if (compute_dipole)
        force->compute_dipole();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else if (method == "caspt2") {

      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      ref = force->ref();
      if (compute_dipole)
        force->compute_dipole();
      out = force->compute(jobtitle, gradinfo);
      force_dipole_ = force->dipole();

    } else {

      numerical_ = true;
      cout << "  There is no analytical gradient available. Numerical gradient will be used." << endl;

    }

    if (export_grad)
      force_export(jobtitle, gradinfo, energyvec, out, export_single);

  }

  if (numerical_) {

    const int target = idata_->get<int>("target", 0);
    const int target2= idata_->get<int>("target2", 1);

    const double dx = idata_->get<double>("dx", 1.0e-3);
    const int nproc = idata_->get<int>("nproc", 1);

    if (jobtitle == "force") {

      auto force = make_shared<FiniteGrad>(input, geom_, ref_, target, dx, nproc);
      out = force->compute();
      ref = force->ref();

    } else if (jobtitle == "nacme") {

      if (method == "casscf") {

        auto force = make_shared<FiniteNacm<CASSCF>>(cinput, geom_, ref_, target, target2, dx, nproc);
        out = force->compute();
        ref = force->ref();

      } else if (method == "caspt2") {

        auto force = make_shared<FiniteNacm<CASPT2Energy>>(cinput, geom_, ref_, target, target2, dx, nproc);
        out = force->compute();
        ref = force->ref();

      }
    }
  }

  ref_ = ref;
  return out;
}


void Force::force_export(const string jobtitle, shared_ptr<const GradInfo> gradinfo, const vector<double> energy, shared_ptr<const GradFile> out, const bool export_single) {
  if (export_single) {
    auto wholemuffle = make_shared<Muffle>("FORCE.out");

    wholemuffle->mute();
    cout << setw(20) << setprecision(10) << energy[gradinfo->target_state()] << endl;
    out->print_export();
    wholemuffle->unmute();
  }

  {
    string mufflename = to_upper(jobtitle) + "_" + to_string(gradinfo->target_state());
    if (jobtitle == "nacme")
      mufflename += ("_" + to_string(gradinfo->target_state2()));
    mufflename += ".out";
    auto gradmuffle = make_shared<Muffle>(mufflename);

    gradmuffle->mute();
    out->print_export();
    gradmuffle->unmute();
  }

  {
    auto enermuffle = make_shared<Muffle>("ENERGY.out");
    enermuffle->mute();
    for (auto i : energy)
      cout << setw(20) << setprecision(10) << i << endl;
    enermuffle->unmute();
  }
}
