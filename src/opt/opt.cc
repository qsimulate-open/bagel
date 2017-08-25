//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: opt.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <functional>
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/opt/optimize.h>
#include <src/opt/opt.h>
#include <src/util/archive.h>
#include <src/grad/hess.h>

using namespace std;
using namespace bagel;

// TODO  Constrained optimization

Opt::Opt(shared_ptr<const PTree> idat, shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : idata_(idat), input_(inp), current_(geom), prev_ref_(ref), iter_(0) {

  auto lastmethod = *idat->get_child("method")->rbegin();
  method_ = to_lower(lastmethod->get<string>("title", ""));

  target_state_ = idat->get<int>("target", 0);
  internal_ = idat->get<bool>("internal", true);
  redundant_ = idat->get<bool>("redundant", false);
  maxiter_ = idat->get<int>("maxiter", 100);
  maxziter_ = idat->get<int>("maxziter", 100);
  scratch_ = idat->get<bool>("scratch", false);
  numerical_ = idat->get<bool>("numerical", false);
  hess_update_ = to_lower(idat->get<string>("hess_update", "flowchart"));
  hess_approx_ = idat->get<bool>("hess_approx", true);
  qmmm_ = idat->get<bool>("qmmm", false);

  constrained_ = idat->get<bool>("constrained", false);
  if (constrained_) {
    if (!internal_ || redundant_) throw runtime_error("Constrained optimization currently only for delocalized internals");
    auto constraints = idat->get_child("constraint");
    for (auto& c : *constraints) {
      constraints_.push_back(make_shared<const OptConstraint>(c));
    }
    cout << "# constraints = " << constraints_.size() << endl;
  }

  explicit_bond_ = idat->get<bool>("explicitbond", false);
  if (explicit_bond_) {
    auto explicit_bonds = idat->get_child("explicit");
    for (auto& e : *explicit_bonds) {
      bonds_.push_back(make_shared<const OptExpBonds>(e));
    }
    cout << endl << "  * Added " << bonds_.size() << " bonds between the non-bonded atoms in overall" << endl;
  }

  if (qmmm_) {
    string qmmm_program = to_lower(idat->get<string>("qmmm_program", "tinker"));
    if (qmmm_program == "tinker") {
      qmmm_driver_ = make_shared<QMMM_Tinker>();
    } else {
      throw runtime_error("QM/MM optimization is only supported with TINKER program");
    }
    // QM/MM does not (and probably, should not) support the internal coordinates
    internal_ = false;
  }

  opttype_ = to_lower(idat->get<string>("opttype", "energy"));
  if (internal_) {
    if (redundant_)
      bmat_red_ = current_->compute_redundant_coordinate();
    else
      bmat_ = current_->compute_internal_coordinate(nullptr, bonds_, constraints_, (opttype_=="transition"));
  }

  // small molecule (atomno < 4) threshold : (1.0e-5, 4.0e-5, 1.0e-6)  (tight in GAUSSIAN and Q-Chem = normal / 30)
  // large molecule              threshold : (3.0e-4, 1.2e-3, 1.0e-6)  (normal in GAUSSIAN and Q-Chem)
  if (current_->natom() < 4 && opttype_ == "energy") {
    thresh_grad_ = idat->get<double>("maxgrad", 0.00001);
    thresh_displ_ = idat->get<double>("maxdisp", 0.00004);
    thresh_echange_ = idat->get<double>("maxchange", 0.000001);
  } else {
    thresh_grad_ = idat->get<double>("maxgrad", 0.0003);
    thresh_displ_ = idat->get<double>("maxdisp", 0.0012);
    thresh_echange_ = idat->get<double>("maxchange", 0.000001);
  }
  maxstep_ = idat->get<double>("maxstep", opttype_ == "energy" ? 0.3 : 0.1);
  algorithm_ = to_lower(idat->get<string>("algorithm", "ef"));
  adaptive_ = idat->get<bool>("adaptive", algorithm_ == "rfo" ? true : false);

  if (opttype_ == "conical" || opttype_ == "meci" || opttype_ == "mdci") {
    // parameters for CI optimizations (Bearpark, Robb, Schlegel)
    target_state2_ = idat->get<int>("target2", 1);
    if (target_state2_ > target_state_) {
      const int tmpstate = target_state_;
      target_state_ = target_state2_;
      target_state2_ = tmpstate;
    }
    nacmtype_ = to_lower(idat->get<string>("nacmtype", "noweight"));
    thielc3_  = idat->get<double>("thielc3", opttype_=="mdci" ? 0.01 : 2.0);
    thielc4_  = idat->get<double>("thielc4", 0.5);
    adaptive_ = false;        // we cannot use it for conical intersection optimization because we do not have a target function
  } else if (opttype_ == "mep") {
    // parameters for MEP calculations (Gonzalez, Schlegel)
    mep_direction_ = idat->get<int>("mep_direction", 1);
    if (hess_approx_) throw runtime_error("MEP calculation should be started with Hessian eigenvectors");
  } else if (opttype_ != "energy" && opttype_ != "transition") {
    throw runtime_error("Optimization type should be: \"energy\", \"transition\", \"conical\" (\"meci\", \"mdci\"), or \"mep\"");
  }
}


void Opt::compute() {
  auto displ = make_shared<XYZFile>(current_->natom());
  size_ = internal_ ? (redundant_? bmat_red_[0]->ndim() : bmat_[0]->mdim()) : current_->natom()*3;

  dispsize_ = max(current_->natom(),int(size_/3+1));
  displ_ = make_shared<XYZFile>(dispsize_);
  grad_ = make_shared<GradFile>(dispsize_);

  auto mep_start = make_shared<XYZFile>(current_->natom());

  if (hess_approx_) {
    cout << "    * Use approximate Hessian for optimization" << endl;
    if (internal_ && !redundant_) hess_ = make_shared<Matrix>(*(bmat_[2]));
    else {
      hess_ = make_shared<Matrix>(size_, size_);
      hess_->unit();
    }
  } else {
    cout << "    * Compute molecular Hessian for optimization" << endl;
    auto hess = make_shared<Hess>(idata_, current_, prev_ref_);
    hess->compute();
    // if internal, we should transform the Hessian according to Eq. (6) in Schlegel
    // dB/dX term currently omitted (reasonable approximation: Handbook of Computational Chemistry, pp. 323--324)
    if (internal_)
      if (redundant_)
        hess_ = make_shared<Matrix>(*bmat_red_[1] % *(hess->hess()) * *bmat_red_[1]);
      else
        hess_ = make_shared<Matrix>(*bmat_[1] % *(hess->hess()) * *bmat_[1]);
    else
      hess_ = hess->hess()->copy();

    copy_n(hess->proj_hess()->element_ptr(0,abs(mep_direction_) - 1), current_->natom()*3, mep_start->data());
  }

  if (opttype_ == "mep") {
    compute_mep(mep_start);
  } else {
    compute_optimize();
  }
}


void Opt::print_header() const {
  if (opttype_ == "energy" || opttype_ == "transition") {
    cout << endl << "  *** Geometry optimization started ***" << endl <<
                              "     iter         energy               grad rms       time"
    << endl << endl;
  } else if (opttype_ == "conical" || opttype_ == "meci") {
    cout << endl << "  *** Conical intersection optimization started ***" << endl <<
                              "     iter         energy             gap energy            grad rms       time"
    << endl << endl;
  } else if (opttype_ == "mdci") {
    cout << endl << "  *** Conical intersection optimization started ***" << endl <<
                              "     iter       distance             gap energy            grad rms       time"
    << endl << endl;
  }
}


void Opt::print_iteration(const double residual, const double param, const double time) const {
  if (opttype_ == "energy" || opttype_ == "transition")
    print_iteration_energy(residual, time);
  else if (opttype_ == "conical" || opttype_ == "meci" || opttype_ == "mdci")
    print_iteration_conical(residual, param, time);
  print_history_molden();
}


void Opt::print_iteration_energy(const double residual, const double time) const {
  cout << setw(7) << iter_ << setw(20) << setprecision(8) << fixed << en_
                           << setw(20) << setprecision(8) << fixed << residual
                           << setw(12) << setprecision(2) << fixed << time << endl;
}


void Opt::print_iteration_conical(const double residual, const double param, const double time) const {
  cout << setw(7) << iter_ << setw(20) << setprecision(8) << fixed << en_
                           << setw(20) << setprecision(8) << fixed << param
                           << setw(20) << setprecision(8) << fixed << residual
                           << setw(12) << setprecision(2) << fixed << time << endl;
}


void Opt::print_history_molden() const {
  const int nopt = prev_en_.size();
  if (nopt != prev_xyz_.size() || nopt != prev_grad_.size())
    throw logic_error("error print_history_molden()");
  stringstream ss;
  ss << " [MOLDEN FORMAT]" << endl;
  ss << " [N_GEO]"         << endl;
  ss << setw(20) << nopt   << endl;
  ss << " [GEOCONV]"       << endl;
  ss << " energy"          << endl;
  for (auto& i : prev_en_)
    ss << scientific << setprecision(20) << i << endl;
  ss << " max-force"       << endl;
  for (auto& i : prev_grad_)
    ss << fixed << setw(15) << setprecision(10) << abs(*max_element(i->begin(), i->end(), [](double x, double y){ return abs(x)<abs(y); })) / au2angstrom__ << endl;
  ss << " rms-force"       << endl;
  for (auto& i : prev_grad_)
    ss << fixed << setw(15) << setprecision(10) << i->rms() / au2angstrom__ << endl;
  ss << " max-step"        << endl;
  for (auto& i : prev_displ_)
    ss << fixed << setw(15) << setprecision(10) << abs(*max_element(i->begin(), i->end(), [](double x, double y){ return abs(x)<abs(y); })) * au2angstrom__ << endl;
  ss << " rms-step"        << endl;
  for (auto& i : prev_displ_)
    ss << fixed << setw(15) << setprecision(10) << i->rms() * au2angstrom__ << endl;

  ss << " [GEOMETRIES] (XYZ)" << endl;
  const int natom = current_->natom();
  for (int i = 0; i != nopt; ++i) {
    ss << setw(4) << natom << endl;
    ss << setw(30) << setprecision(20) << prev_en_.at(i) << endl;
    for (int j = 0; j != natom; ++j) {
      string name = current_->atoms(j)->name();
      name[0] = toupper(name[0]);
      ss << name << setw(20) << setprecision(10) << prev_xyz_.at(i)->element(0, j) * au2angstrom__
                 << setw(20) << setprecision(10) << prev_xyz_.at(i)->element(1, j) * au2angstrom__
                 << setw(20) << setprecision(10) << prev_xyz_.at(i)->element(2, j) * au2angstrom__ << endl;
    }
  }

  ss << " [FORCES]" << endl;
  for (int i = 0; i != nopt; ++i) {
    ss << "point" << setw(4) << i << endl;
    ss << setw(4) << natom << endl;
    for (int j = 0; j != natom; ++j) {
      ss << setw(20) << setprecision(10) << prev_grad_.at(i)->element(0, j) / au2angstrom__
         << setw(20) << setprecision(10) << prev_grad_.at(i)->element(1, j) / au2angstrom__
         << setw(20) << setprecision(10) << prev_grad_.at(i)->element(2, j) / au2angstrom__ << endl;
    }
  }

  ofstream fs("opt_history.molden");
  fs << ss.str();
}

