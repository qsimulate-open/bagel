//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradinfo.h
// Copyright (C) 2017 Toru Shiozaki
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


#ifndef __SRC_GRAD_GRADINFO_H
#define __SRC_GRAD_GRADINFO_H

#include <src/grad/nacmtype.h>

namespace bagel {

class GradInfo {
  protected:
    int target_state_;
    int target_state2_;
    int maxziter_;
    bool density_print_;
    bool cider_eval_;
    std::shared_ptr<NacmType> nacmtype_;
    std::shared_ptr<PTree> moprint_info_;

  public:
    GradInfo() : target_state_(0), target_state2_(1), maxziter_(100), density_print_(false), cider_eval_(true), nacmtype_(std::make_shared<NacmType>("full")), moprint_info_(std::make_shared<PTree>()) { }
    GradInfo(std::shared_ptr<const GradInfo> o, const int target) {
      target_state_ = target;
      target_state2_ = -1;
      maxziter_ = o->maxziter();
      density_print_ = o->density_print();
      cider_eval_ = o->cider_eval();
      nacmtype_ = o->nacmtype();
      moprint_info_ = o->moprint_info();
    }
    GradInfo(std::shared_ptr<const PTree> idat, const bool opt = false) {
      nacmtype_ = std::make_shared<NacmType>(to_lower(idat->get<std::string>("nacmtype", opt ? "noweight" : "full")));
      target_state_ = idat->get<int>("target", 0);
      target_state2_ = idat->get<int>("target2", 1);
      maxziter_ = idat->get<int>("maxziter", 100);
      density_print_ = idat->get<bool>("density_print", false);
      cider_eval_ = idat->get<bool>("ciderivative", true);
      moprint_info_ = idat->get_child_optional("moprint");
      if (!moprint_info_) moprint_info_ = std::make_shared<PTree>();
    }

    int target_state() const { return target_state_; }
    int target_state2() const { return target_state2_; }
    int maxziter() const { return maxziter_; }
    bool density_print() const { return density_print_; }
    bool cider_eval() const { return cider_eval_; }
    std::shared_ptr<NacmType> nacmtype() const { return nacmtype_; }
    std::shared_ptr<PTree> moprint_info() const { return moprint_info_; }

};

}

#endif
