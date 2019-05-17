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
  : idata_(idat), input_(inp), current_(geom), prev_ref_(ref) {

  auto lastmethod = *idat->get_child("method")->rbegin();
  method_ = to_lower(lastmethod->get<string>("title", ""));

  optinfo_ = make_shared<const OptInfo>(idat, geom);

  if (optinfo_->qmmm()) {
    string qmmm_program = to_lower(idat->get<string>("qmmm_program", "tinker"));
    if (qmmm_program == "tinker") {
      qmmm_driver_ = make_shared<const QMMM_Tinker>();
    } else {
      throw runtime_error("QM/MM optimization is only supported with TINKER program");
    }
  }

  if (optinfo_->internal()) {
    if (optinfo_->redundant())
      tie(bondlist_, bmat_red_) = current_->compute_redundant_coordinate(bondlist_);
    else
      bmat_ = current_->compute_internal_coordinate(nullptr, optinfo_->bonds(), optinfo_->opttype()->is_transition());
  }

  maxstep_ = idat->get<double>("maxstep", optinfo_->opttype()->is_energy() || optinfo_->opttype()->is_transition() ? 0.3 : 0.1);

}


void Opt::compute() {
  auto displ = make_shared<XYZFile>(current_->natom());
  size_ = optinfo()->internal() ? (optinfo()->redundant()? bmat_red_[0]->ndim() : bmat_[0]->mdim()) : current_->natom()*3;

  dispsize_ = max(current_->natom(), int(size_ / 3 + 1));
  displ_ = make_shared<XYZFile>(dispsize_);
  grad_ = make_shared<GradFile>(dispsize_);

  auto mep_start = make_shared<XYZFile>(current_->natom());

  if (optinfo()->hess_approx()) {
    cout << "    * Use approximate Hessian for optimization" << endl;
    if (optinfo()->internal()) {
      if (optinfo()->redundant()) {
        hess_ = make_shared<Matrix>(*(bmat_red_[3]));
      } else {
        hess_ = make_shared<Matrix>(*(bmat_[2]));
      }
    } else {
      hess_ = make_shared<Matrix>(size_, size_);
      hess_->unit();
    }
  } else {
    cout << "    * Compute molecular Hessian for optimization" << endl;
    auto hess = make_shared<Hess>(idata_, current_, prev_ref_);
    hess->compute();
    // if internal, we should transform the Hessian according to Eq. (6) in Schlegel
    // dB/dX term currently omitted (reasonable approximation: Handbook of Computational Chemistry, pp. 323--324)
    if (optinfo()->internal()) {
      if (optinfo()->redundant())
        hess_ = make_shared<Matrix>(*bmat_red_[1] % *(hess->hess()) * *bmat_red_[1]);
      else
        hess_ = make_shared<Matrix>(*bmat_[1] % *(hess->hess()) * *bmat_[1]);
    } else {
      hess_ = hess->hess()->copy();
    }

    if (optinfo()->opttype()->is_mep()) {
      if (optinfo()->mep_direction() == 0) {
        mep_start->zero();
      } else {
        copy_n(hess->proj_hess()->element_ptr(0, abs(optinfo()->mep_direction()) - 1), current_->natom() * 3, mep_start->data());
      }
    }
  }

  if (optinfo_->opttype()->is_mep()) {
    compute_mep(mep_start);
  } else {
    compute_optimize();
  }
}


void Opt::print_header() const {
  if (optinfo()->opttype()->is_energy() || optinfo()->opttype()->is_transition()) {
    cout << endl << "  *** Geometry optimization started ***" << endl <<
                              "     iter         energy               grad rms       time"
    << endl << endl;
  } else if (optinfo()->opttype()->is_mdci()) {
    cout << endl << "  *** Conical intersection optimization started ***" << endl <<
                              "     iter       distance             gap energy            grad rms       time"
    << endl << endl;
  } else if (optinfo()->opttype()->is_conical()) {
    cout << endl << "  *** Conical intersection optimization started ***" << endl <<
                              "     iter         energy             gap energy            grad rms       time"
    << endl << endl;
  }
}


void Opt::print_iteration(const int iter, const double residual, const double param, const double time) const {
  if (optinfo()->opttype()->is_conical()) {
    print_iteration_conical(iter, residual, param, time);
  } else {
    print_iteration_energy(iter, residual, time);
  }
  print_history_molden();

  if (optinfo_->molden()) {
    stringstream ss; ss << "geom" << iter << ".molden";
    MoldenOut m(ss.str());
    m << prev_ref_->geom();
    m << prev_ref_;
  }
}


void Opt::print_iteration_energy(const int iter, const double residual, const double time) const {
  cout << setw(7) << iter << setw(20) << setprecision(8) << fixed << en_
                          << setw(20) << setprecision(8) << fixed << residual
                          << setw(12) << setprecision(2) << fixed << time << endl;
}


void Opt::print_iteration_conical(const int iter, const double residual, const double param, const double time) const {
  cout << setw(7) << iter << setw(20) << setprecision(8) << fixed << en_
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

