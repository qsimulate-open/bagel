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
  const bool export_grad = idata_->get<bool>("export", false);
  const bool export_single = idata_->get<bool>("export_single", false);

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
  vector<double> energyvec;

  if (jobtitle == "forces") {

    // list of gradients (and NACMEs) to be evaluated. maxziter, nacmtype, target, target2 can be specified
    auto joblist = idata_->get_child("grads");

    if (method == "casscf") {

      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();

      for (auto& m : *joblist) {
        const string mtitle= to_lower(m->get<string>("title", "force"));
        const int target  = m->get<int>("target", 0);
        const int target2 = m->get<int>("target2", 1);
        const int maxziter = m->get<int>("maxziter", 100);
        const int nacmtype = m->get<int>("nacmtype", 0);
        out = force->compute(mtitle, target, maxziter, target2, nacmtype);

        if (export_grad)
          force_export(export_single, target, target2, energyvec, mtitle, out);
      }

    } else if (method == "caspt2") {

      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_);
      energyvec = force->energyvec();

      for (auto& m : *joblist) {
        const string mtitle= to_lower(m->get<string>("title", "force"));
        const int target  = m->get<int>("target", 0);
        const int target2 = m->get<int>("target2", 1);
        const int maxziter = m->get<int>("maxziter", 100);
        const int nacmtype = m->get<int>("nacmtype", 0);
        out = force->compute(mtitle, target, maxziter, target2, nacmtype);

        if (export_grad)
          force_export(export_single, target, target2, energyvec, mtitle, out);
      }

    } else {

      throw runtime_error("multiple force calculations can be only done with CASSCF and CASPT2");

    }

  } else if (!numerical_) {

    const int target = idata_->get<int>("target", 0);
    const int target2= idata_->get<int>("target2", 1);

    if (method == "uhf") {

      auto force = make_shared<GradEval<UHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "rohf") {

      auto force = make_shared<GradEval<ROHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "hf") {

      auto force = make_shared<GradEval<RHF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "ks") {

      auto force = make_shared<GradEval<KS>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "dhf") {

      auto force = make_shared<GradEval<Dirac>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target);
      ref = force->ref();

    } else if (method == "mp2") {

      const int maxziter = idata_->get<int>("maxziter", 100);
      auto force = make_shared<GradEval<MP2Grad>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target, maxziter);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "casscf") {

      const int maxziter = idata_->get<int>("maxziter", 100);
      const int nacmtype = idata_->get<int>("nacmtype", 0);
      auto force = make_shared<GradEval<CASSCF>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target, maxziter, target2, nacmtype);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else if (method == "caspt2") {

      const int maxziter = idata_->get<int>("maxziter", 100);
      const int nacmtype = idata_->get<int>("nacmtype", 0);
      auto force = make_shared<GradEval<CASPT2Grad>>(cinput, geom_, ref_);
      energyvec = force->energyvec();
      out = force->compute(jobtitle, target, maxziter, target2, nacmtype);
      ref = force->ref();
      force_dipole_ = force->dipole();

    } else {

        numerical_ = true;
        cout << "  It seems like no analytical gradient method available; moving to finite difference " << endl;

    }

    if (export_grad)
      force_export(export_single, target, target2, energyvec, jobtitle, out);

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

  return out;
}


void Force::force_export(const bool export_single, const int target, const int target2, const vector<double> energy, const string jobtitle, const shared_ptr<GradFile> out) {
  if (export_single) {
    shared_ptr<Muffle> wholemuffle;
    wholemuffle = make_shared<Muffle>("FORCE.out");
    wholemuffle->mute();
    cout << setw(20) << setprecision(10) << energy[target] << endl;
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
  for (auto i : energy)
    cout << setw(20) << setprecision(10) << i << endl;
  enermuffle->unmute();

}
