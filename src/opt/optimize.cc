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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

Optimize::Optimize(const boost::property_tree::ptree& idata, std::shared_ptr<const Geometry> g) : idata_(idata), geom_(g) {
  maxiter_ = idata.get<int>("maxiter", 100);

}


void Optimize::compute() {
  const boost::property_tree::ptree methodblock = idata_.get_child("method");
  string method = methodblock.get<string>("title", ""); 
  if (method.empty()) throw std::runtime_error("title is missing in one of the input blocks (opt)");

  if (method == "df-uhf" || method == "uhf") {

    std::shared_ptr<Opt<UHF>> opt(new Opt<UHF>(idata_, methodblock, geom_));
    for (int i = 0; i != maxiter_; ++i)
      if (opt->next()) break;
    geom_ = opt->geometry();

  } else if (method == "df-rohf" || method == "rohf") {

    std::shared_ptr<Opt<ROHF>> opt(new Opt<ROHF>(idata_, methodblock, geom_));
    for (int i = 0; i != maxiter_; ++i)
      if (opt->next()) break;
    geom_ = opt->geometry();

  } else if (method == "df-hf") {

    std::shared_ptr<Opt<SCF<1>>> opt(new Opt<SCF<1>>(idata_, methodblock, geom_));
    for (int i = 0; i != maxiter_; ++i)
      if (opt->next()) break;
    geom_ = opt->geometry();

  } else if (method == "df-ks") {

    std::shared_ptr<Opt<KS>> opt(new Opt<KS>(idata_, methodblock, geom_));
    for (int i = 0; i != maxiter_; ++i)
      if (opt->next()) break;
    geom_ = opt->geometry();

  } else if (method == "mp2") {

    std::shared_ptr<Opt<MP2Grad>> opt(new Opt<MP2Grad>(idata_, methodblock, geom_));
    for (int i = 0; i != maxiter_; ++i)
      if (opt->next()) break;
    geom_ = opt->geometry();

  } else if (method == "casscf") {
    std::string algorithm = methodblock.get<std::string>("algorithm", "");
    // in case of SS-CASSCF
    if (methodblock.get<int>("nstate", 1) == 1) {
      if (algorithm == "superci" || algorithm == "") {
        std::shared_ptr<Opt<SuperCI>> opt(new Opt<SuperCI>(idata_, methodblock, geom_));
        for (int i = 0; i != maxiter_; ++i)
          if (opt->next()) break;
        geom_ = opt->geometry();
      } else if (algorithm == "werner" || algorithm == "knowles") {
        std::shared_ptr<Opt<WernerKnowles>> opt(new Opt<WernerKnowles>(idata_, methodblock, geom_));
        for (int i = 0; i != maxiter_; ++i)
          if (opt->next()) break;
        geom_ = opt->geometry();
      } else {
        throw std::runtime_error("unknown CASSCF algorithm specified.");
      }
    // in case of SA-CASSCF
    } else {
      if (algorithm == "superci" || algorithm == "") {
        std::shared_ptr<Opt<SuperCIGrad>> opt(new Opt<SuperCIGrad>(idata_, methodblock, geom_));
        for (int i = 0; i != maxiter_; ++i)
          if (opt->next()) break;
        geom_ = opt->geometry();
      } else {
        throw std::runtime_error("unknown CASSCF algorithm specified.");
      }
    }
  }

}
