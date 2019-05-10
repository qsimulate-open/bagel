//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: optinfo.h
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


#ifndef __SRC_OPT_OPTINFO_H
#define __SRC_OPT_OPTINFO_H

#include <src/grad/gradinfo.h>
#include <src/opt/constraint.h>

namespace bagel {

enum OptTargetType { energy, meci, mdci, transition, mep };
enum OptAlgorithmType { ef, rfo, nr };
enum HessUpdateType { flowchart, bfgs, psb, sr1, bofill, noupdate };


class HessUpdate {
  protected:
    HessUpdateType type_;

  public:
    HessUpdate() : type_(HessUpdateType::flowchart) { }
    HessUpdate(std::string input_hess_update) {
      if (input_hess_update == "flowchart") {
        type_ = HessUpdateType::flowchart;
      } else if (input_hess_update == "bfgs") {
        type_ = HessUpdateType::bfgs;
      } else if (input_hess_update == "psb") {
        type_ = HessUpdateType::psb;
      } else if (input_hess_update == "sr1") {
        type_ = HessUpdateType::sr1;
      } else if (input_hess_update == "bofill") {
        type_ = HessUpdateType::bofill;
      } else if (input_hess_update == "noupdate") {
        type_ = HessUpdateType::noupdate;
      } else {
        throw std::logic_error ("Available hess_update: \"flowchart\", \"bfgs\", \"psb\", \"sr1\", \"bofill\", \"noupdate\".");
      }
    }

    bool is_flowchart() const { return type_ == HessUpdateType::flowchart; }
    bool is_bfgs() const { return type_ == HessUpdateType::bfgs; }
    bool is_psb() const { return type_ == HessUpdateType::psb; }
    bool is_sr1() const { return type_ == HessUpdateType::sr1; }
    bool is_bofill() const { return type_ == HessUpdateType::bofill; }
    bool is_noupdate() const { return type_ == HessUpdateType::noupdate; }
};


class OptAlgorithms {
  protected:
    OptAlgorithmType type_;

  public:
    OptAlgorithms() : type_(OptAlgorithmType::ef) { }
    OptAlgorithms(std::shared_ptr<const PTree> idat) {
      const std::string input_algorithm = to_lower(idat->get<std::string>("algorithm", "ef"));
      if (input_algorithm == "ef") {
        type_ = OptAlgorithmType::ef;
      } else if (input_algorithm == "rfo") {
        type_ = OptAlgorithmType::rfo;
      } else if (input_algorithm == "nr") {
        type_ = OptAlgorithmType::nr;
      } else {
        throw std::logic_error ("Available algorithms: \"ef\", \"rfo\", or \"nr\".");
      }
    }

    bool is_ef() const { return type_ == OptAlgorithmType::ef; }
    bool is_rfo() const { return type_ == OptAlgorithmType::rfo; }
    bool is_nr() const { return type_ == OptAlgorithmType::nr; }
};


class OptType {
  protected:
    OptTargetType type_;
    bool conical_;

  public:
    OptType() : type_(OptTargetType::energy), conical_(false) { }
    OptType(std::string input_opttype) {
      if (input_opttype == "energy") {
        type_ = OptTargetType::energy;
        conical_ = false;
      } else if (input_opttype == "meci" || input_opttype == "conical") {
        type_ = OptTargetType::meci;
        conical_ = true;
      } else if (input_opttype == "mdci") {
        type_ = OptTargetType::mdci;
        conical_ = true;
      } else if (input_opttype == "transition") {
        type_ = OptTargetType::transition;
        conical_ = false;
      } else if (input_opttype == "mep") {
        type_ = OptTargetType::mep;
        conical_ = false;
      } else {
        throw std::logic_error ("Available opttypes: \"energy\", \"conical\", \"meci\", \"mdci\", \"transition\", or \"mep\".");
      }
    }

    bool is_conical() const { return conical_; }

    bool is_energy() const { return type_ == OptTargetType::energy; }
    bool is_meci() const { return type_ == OptTargetType::meci; }
    bool is_mdci() const { return type_ == OptTargetType::mdci; }
    bool is_transition() const { return type_ == OptTargetType::transition; }
    bool is_mep() const { return type_ == OptTargetType::mep; }
};


class OptInfo : public GradInfo {
  protected:
    std::shared_ptr<OptType> opttype_;
    std::shared_ptr<OptAlgorithms> algorithm_;
    std::shared_ptr<HessUpdate> hessupdate_;

    bool qmmm_;

    int maxiter_;
    double thresh_grad_;
    double thresh_displ_;
    double thresh_echange_;

    bool scratch_;
    bool numerical_;

    bool internal_;
    bool redundant_;

    bool adaptive_;
    bool hess_approx_;
    double thielc3_;
    double thielc4_;
    int hess_recalc_freq_;

    int mep_direction_;

    bool explicit_bond_;
    std::vector<std::shared_ptr<const OptExpBonds>> bonds_;

    bool molden_;

  public:
    OptInfo(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom) : GradInfo(idat, /*opt=*/true) {
      opttype_ = std::make_shared<OptType>(to_lower(idat->get<std::string>("opttype", "energy")));
      algorithm_ = std::make_shared<OptAlgorithms>(idat);
      hessupdate_ = std::make_shared<HessUpdate>(to_lower(idat->get<std::string>("hess_update", opttype_->is_transition() ? "bofill" : "flowchart")));

      internal_ = idat->get<bool>("internal", true);
      redundant_ = idat->get<bool>("redundant", false);
      maxiter_ = idat->get<int>("maxiter", 100);
      scratch_ = idat->get<bool>("scratch", false);
      numerical_ = idat->get<bool>("numerical", false);
      hess_approx_ = idat->get<bool>("hess_approx", opttype_->is_mep() ? false : true);
      hess_recalc_freq_ = idat->get<int>("hess_recalc_freq", 5);

      explicit_bond_ = idat->get<bool>("explicitbond", false);
      if (explicit_bond_) {
        auto explicit_bonds = idat->get_child("explicit");
        for (auto& e : *explicit_bonds) {
          bonds_.push_back(std::make_shared<const OptExpBonds>(e));
        }
        std::cout << std::endl << "  * Added " << bonds_.size() << " bonds between the non-bonded atoms" << std::endl;
      }

      molden_ = idat->get<bool>("molden", false);

      qmmm_ = idat->get<bool>("qmmm", false);
      if (qmmm_)
        internal_ = false;

      // small molecule (atomno < 4) threshold : (1.0e-5, 4.0e-5, 1.0e-6)  (tight in GAUSSIAN and Q-Chem = normal / 30)
      // large molecule              threshold : (3.0e-4, 1.2e-3, 1.0e-6)  (normal in GAUSSIAN and Q-Chem)
      if (geom->natom() < 4 && opttype_->is_energy()) {
        thresh_grad_ = idat->get<double>("maxgrad", 0.00001);
        thresh_displ_ = idat->get<double>("maxdisp", 0.00004);
        thresh_echange_ = idat->get<double>("maxchange", 0.000001);
      } else {
        thresh_grad_ = idat->get<double>("maxgrad", 0.0003);
        thresh_displ_ = idat->get<double>("maxdisp", 0.0012);
        thresh_echange_ = idat->get<double>("maxchange", 0.000001);
      }

      adaptive_ = idat->get<bool>("adaptive", algorithm_->is_rfo() ? true : false);

      if (opttype_->is_conical()) {
        // parameters for CI optimizations (Bearpark, Robb, Schlegel)
        if (target_state2_ > target_state_) {
          const int tmpstate = target_state_;
          target_state_ = target_state2_;
          target_state2_ = tmpstate;
        }
        nacmtype_ = std::make_shared<NacmType>(to_lower(idat->get<std::string>("nacmtype", "noweight")));
        thielc3_  = idat->get<double>("thielc3", opttype_->is_mdci() ? 0.01 : 2.0);
        thielc4_  = idat->get<double>("thielc4", 0.5);
        adaptive_ = false;        // we cannot use it for conical intersection optimization because we do not have a target function
      } else {
        // initialize the values
        nacmtype_ = std::make_shared<NacmType>();
        thielc3_ = 2.0;
        thielc4_ = 0.5;
      }

      if (opttype_->is_mep()) {
        // parameters for MEP calculations (Gonzalez, Schlegel)
        mep_direction_ = idat->get<int>("mep_direction", 1);
        if (hess_approx_ && mep_direction_ != 0)
          throw std::runtime_error("MEP calculation should be started with Hessian eigenvectors");
      } else {
        // initialize the values
        mep_direction_ = 0;
      }
    }


    std::shared_ptr<OptType> opttype() const { return opttype_; }
    std::shared_ptr<OptAlgorithms> algorithm() const { return algorithm_; }
    std::shared_ptr<HessUpdate> hessupdate() const { return hessupdate_; }

    bool qmmm() const { return qmmm_; }

    int maxiter() const { return maxiter_; }
    double thresh_grad() const { return thresh_grad_; }
    double thresh_displ() const { return thresh_displ_; }
    double thresh_echange() const { return thresh_echange_; }

    bool scratch() const { return scratch_; }
    bool numerical() const { return numerical_; }

    bool internal() const { return internal_; }
    bool redundant() const { return redundant_; }

    bool adaptive() const { return adaptive_; }
    bool hess_approx() const { return hess_approx_; }
    double thielc3() const { return thielc3_; }
    double thielc4() const { return thielc4_; }
    int hess_recalc_freq() const { return hess_recalc_freq_; }

    int mep_direction() const { return mep_direction_; }

    bool explicit_bond() const { return explicit_bond_; }
    std::vector<std::shared_ptr<const OptExpBonds>> bonds() const { return bonds_; }

    bool molden() const { return molden_; }
};

}

#endif
