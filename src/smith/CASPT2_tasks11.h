//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks11.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#ifndef __SRC_SMITH_CASPT2_TASKS11_H
#define __SRC_SMITH_CASPT2_TASKS11_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task500 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task500(std::shared_ptr<TATensor<double,2>> I687, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I687), ta1_(Gamma60), ta2_(t2) { }
};

class Task501 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task501(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I689)
   : ta0_(proj), ta1_(I689) { }
};

class Task502 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task502(std::shared_ptr<TATensor<double,2>> I689, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,2>> I690)
   : ta0_(I689), ta1_(Gamma38), ta2_(I690) { }
};

class Task503 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task503(std::shared_ptr<TATensor<double,2>> I690, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I690), ta1_(t2) { }
};

class Task504 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task504(std::shared_ptr<TATensor<double,2>> I690, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I690), ta1_(t2) { }
};

class Task505 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task505(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I701)
   : ta0_(proj), ta1_(I701) { }
};

class Task506 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task506(std::shared_ptr<TATensor<double,2>> I701, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I702)
   : ta0_(I701), ta1_(t2), ta2_(I702) { }
};

class Task507 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task507(std::shared_ptr<TATensor<double,4>> I702, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I702), ta1_(Gamma38), ta2_(t2) { }
};

class Task508 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task508(std::shared_ptr<TATensor<double,2>> I701, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I705)
   : ta0_(I701), ta1_(t2), ta2_(I705) { }
};

class Task509 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task509(std::shared_ptr<TATensor<double,4>> I705, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I705), ta1_(Gamma38), ta2_(t2) { }
};

class Task510 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task510(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I707)
   : ta0_(proj), ta1_(I707) { }
};

class Task511 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task511(std::shared_ptr<TATensor<double,2>> I707, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I708)
   : ta0_(I707), ta1_(t2), ta2_(I708) { }
};

class Task512 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task512(std::shared_ptr<TATensor<double,4>> I708, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I708), ta1_(Gamma38), ta2_(t2) { }
};

class Task513 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task513(std::shared_ptr<TATensor<double,2>> I707, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I711)
   : ta0_(I707), ta1_(t2), ta2_(I711) { }
};

class Task514 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task514(std::shared_ptr<TATensor<double,4>> I711, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I711), ta1_(Gamma38), ta2_(t2) { }
};

class Task515 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task515(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I713)
   : ta0_(proj), ta1_(I713) { }
};

class Task516 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task516(std::shared_ptr<TATensor<double,2>> I713, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I714)
   : ta0_(I713), ta1_(t2), ta2_(I714) { }
};

class Task517 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task517(std::shared_ptr<TATensor<double,4>> I714, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I714), ta1_(Gamma38), ta2_(t2) { }
};

class Task518 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task518(std::shared_ptr<TATensor<double,2>> I713, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I717)
   : ta0_(I713), ta1_(t2), ta2_(I717) { }
};

class Task519 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task519(std::shared_ptr<TATensor<double,4>> I717, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I717), ta1_(Gamma38), ta2_(t2) { }
};

class Task520 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task520(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I719)
   : ta0_(proj), ta1_(I719) { }
};

class Task521 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task521(std::shared_ptr<TATensor<double,2>> I719, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I720)
   : ta0_(I719), ta1_(t2), ta2_(I720) { }
};

class Task522 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task522(std::shared_ptr<TATensor<double,4>> I720, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I720), ta1_(Gamma60), ta2_(t2) { }
};

class Task523 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task523(std::shared_ptr<TATensor<double,2>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task524 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task524(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I734)
   : ta0_(proj), ta1_(I734) { }
};

class Task525 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task525(std::shared_ptr<TATensor<double,2>> I734, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I734), ta1_(Gamma12), ta2_(t2) { }
};

class Task526 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task526(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I736)
   : ta0_(proj), ta1_(I736) { }
};

class Task527 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task527(std::shared_ptr<TATensor<double,2>> I736, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I737)
   : ta0_(I736), ta1_(Gamma38), ta2_(I737) { }
};

class Task528 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task528(std::shared_ptr<TATensor<double,4>> I737, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I737), ta1_(t2) { }
};

class Task529 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task529(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I740)
   : ta0_(proj), ta1_(I740) { }
};

class Task530 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task530(std::shared_ptr<TATensor<double,2>> I740, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I740), ta1_(Gamma60), ta2_(t2) { }
};

class Task531 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task531(std::shared_ptr<TATensor<double,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task532 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task532(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I742)
   : ta0_(proj), ta1_(I742) { }
};

class Task533 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task533(std::shared_ptr<TATensor<double,4>> I742, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I742), ta1_(Gamma92), ta2_(t2) { }
};

class Task534 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task534(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I744)
   : ta0_(proj), ta1_(I744) { }
};

class Task535 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task535(std::shared_ptr<TATensor<double,4>> I744, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I744), ta1_(Gamma6), ta2_(t2) { }
};

class Task536 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task536(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I746)
   : ta0_(proj), ta1_(I746) { }
};

class Task537 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task537(std::shared_ptr<TATensor<double,4>> I746, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> I747)
   : ta0_(I746), ta1_(Gamma16), ta2_(I747) { }
};

class Task538 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task538(std::shared_ptr<TATensor<double,4>> I747, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I747), ta1_(t2) { }
};

class Task539 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task539(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I750)
   : ta0_(proj), ta1_(I750) { }
};

class Task540 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task540(std::shared_ptr<TATensor<double,4>> I750, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I750), ta1_(Gamma32), ta2_(t2) { }
};

class Task541 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task541(std::shared_ptr<TATensor<double,4>> I750, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I750), ta1_(Gamma35), ta2_(t2) { }
};

class Task542 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task542(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I754)
   : ta0_(proj), ta1_(I754) { }
};

class Task543 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task543(std::shared_ptr<TATensor<double,4>> I754, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> I755)
   : ta0_(I754), ta1_(Gamma35), ta2_(I755) { }
};

class Task544 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task544(std::shared_ptr<TATensor<double,4>> I755, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I755), ta1_(t2) { }
};

class Task545 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task545(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I758)
   : ta0_(proj), ta1_(I758) { }
};

class Task546 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task546(std::shared_ptr<TATensor<double,4>> I758, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I758), ta1_(Gamma59), ta2_(t2) { }
};

class Task547 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task547(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I760)
   : ta0_(proj), ta1_(I760) { }
};

class Task548 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task548(std::shared_ptr<TATensor<double,4>> I760, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I760), ta1_(t2) { }
};

class Task549 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task549(std::shared_ptr<TATensor<double,4>> proj, std::shared_ptr<TATensor<double,4>> I762)
   : ta0_(proj), ta1_(I762) { }
};


}
}
}
#endif
#endif

