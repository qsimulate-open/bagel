//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: get_grad.cc
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
#include <src/wfn/construct_method.h>
#include <src/opt/opt.h>

using namespace std;
using namespace bagel;

shared_ptr<GradFile> Opt::get_cigrad_bearpark(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) {

  // Find conical intersection by gradient projection method (Schlegel)

  auto out = make_shared<GradFile>(current_->natom());
  shared_ptr<GradFile> cgrad1;
  shared_ptr<GradFile> cgrad2;
  shared_ptr<GradFile> x2;
  double en1 = 0.0, en2 = 0.0;

  if (method_ == "casscf") {
    GradEval<CASSCF> eval1(cinput, current_, ref, target_state_);
    cgrad1 = eval1.compute();
    shared_ptr<const Reference> refs = eval1.ref();
    en2 = eval1.energy();

    GradEval<CASSCF> eval2(cinput, current_, refs, target_state2_);
    cgrad2 = eval2.compute();
    refs = eval1.ref();
    en1 = eval2.energy();

    NacmEval<CASSCF> evaln(cinput, current_ ,refs, target_state2_, target_state_, nacmtype_);
    x2 = evaln.compute();
  } else if (method_ == "caspt2") {
    GradEval<CASPT2Grad> eval1(cinput, current_, ref, target_state_);
    cgrad1 = eval1.compute();
    shared_ptr<const Reference> refs = eval1.ref();
    en2 = eval1.energy();

    GradEval<CASPT2Grad> eval2(cinput, current_, refs, target_state2_);
    cgrad2 = eval2.compute();
    refs = eval1.ref();
    en1 = eval2.energy();

    NacmEval<CASPT2Nacm> evaln(cinput, current_, refs, target_state2_, target_state_, nacmtype_);
    x2 = evaln.compute();
  } else {
    throw logic_error ("Conical intersection search currently only available for CASSCF or CASPT2"); 
  }

  auto x1 = std::make_shared<GradFile>(*cgrad2 - *cgrad1);
  auto xf = std::make_shared<GradFile>(*x1);
  auto xg = std::make_shared<GradFile>(*cgrad1);
  const double en  = en2 - en1;
  double x1norm = x1->norm();
  double x2norm = x2->norm();
  xf->scale(-2.0 * en / x1norm);
  x1->scale(1.0 / x1norm);
  x2->scale(1.0 / x2norm);

  auto proj = make_shared<Matrix>(size_, size_);
  proj->unit();
  auto ppt = make_shared<Matrix>(size_,size_);
  auto qqt = make_shared<Matrix>(size_,size_);
  dger_(size_,size_,-1.0,x1->data(),1,x1->data(),1,ppt->data(),size_);
  dger_(size_,size_,-1.0,x2->data(),1,x2->data(),1,qqt->data(),size_);
  *proj = *proj + *ppt + *qqt;
  xg->transform(proj, /*transpose=*/false);

  *out = thielc3_ * (*xf * thielc4_ + *xg * (1.0 - thielc4_));
  en_ = en2;
  egap_ = en;

  return out;
}

shared_ptr<GradFile> Opt::get_grad_energy(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) {
  auto out = make_shared<GradFile>(current_->natom());

  if (method_ == "uhf") {

    GradEval<UHF> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "rohf") {

    GradEval<ROHF> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "hf") {

    GradEval<RHF> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "ks") {

    GradEval<KS> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "dhf") {

    GradEval<Dirac> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "mp2") {

    GradEval<MP2Grad> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "casscf") {

    GradEval<CASSCF> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  } else if (method_ == "caspt2") {

    GradEval<CASPT2Grad> eval(cinput, current_, ref, target_state_);
    out = eval.compute();
    prev_ref_ = eval.ref();
    en_ = eval.energy();
    
  }
  return out;
}

shared_ptr<GradFile> Opt::get_grad(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) {
  auto out = make_shared<GradFile>(current_->natom());
  
  if (opttype_ == "conical") out = get_cigrad_bearpark(cinput, ref);
  else out = get_grad_energy(cinput, ref);

  return out;
}

