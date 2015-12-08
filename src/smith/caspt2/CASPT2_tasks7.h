//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks7.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS7_H
#define __SRC_SMITH_CASPT2_TASKS7_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task300 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task300(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I444)
   : ta0_(proj), ta1_(I444) { }
};

class Task301 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task301(std::shared_ptr<TATensor<double,2>> I444, std::shared_ptr<TATensor<double,6>> Gamma160, std::shared_ptr<TATensor<double,4>> I445)
   : ta0_(I444), ta1_(Gamma160), ta2_(I445) { }
};

class Task302 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task302(std::shared_ptr<TATensor<double,4>> I445, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I445), ta1_(t2) { }
};

class Task303 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task303(std::shared_ptr<TATensor<double,2>> I444, std::shared_ptr<TATensor<double,6>> Gamma191, std::shared_ptr<TATensor<double,4>> I538)
   : ta0_(I444), ta1_(Gamma191), ta2_(I538) { }
};

class Task304 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task304(std::shared_ptr<TATensor<double,4>> I538, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I538), ta1_(t2) { }
};

class Task305 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task305(std::shared_ptr<TATensor<double,2>> I444, std::shared_ptr<TATensor<double,6>> Gamma194, std::shared_ptr<TATensor<double,4>> I547)
   : ta0_(I444), ta1_(Gamma194), ta2_(I547) { }
};

class Task306 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task306(std::shared_ptr<TATensor<double,4>> I547, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I547), ta1_(t2) { }
};

class Task307 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task307(std::shared_ptr<TATensor<double,4>> I547, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I547), ta1_(t2) { }
};

class Task308 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task308(std::shared_ptr<TATensor<double,4>> I547, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I547), ta1_(t2) { }
};

class Task309 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task309(std::shared_ptr<TATensor<double,2>> I444, std::shared_ptr<TATensor<double,6>> Gamma252, std::shared_ptr<TATensor<double,4>> I729)
   : ta0_(I444), ta1_(Gamma252), ta2_(I729) { }
};

class Task310 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task310(std::shared_ptr<TATensor<double,4>> I729, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I729), ta1_(t2) { }
};

class Task311 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task311(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I447)
   : ta0_(proj), ta1_(I447) { }
};

class Task312 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task312(std::shared_ptr<TATensor<double,2>> I447, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I448)
   : ta0_(I447), ta1_(t2), ta2_(I448) { }
};

class Task313 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task313(std::shared_ptr<TATensor<double,4>> I448, std::shared_ptr<TATensor<double,4>> Gamma92, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I448), ta1_(Gamma92), ta2_(t2) { }
};

class Task314 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task314(std::shared_ptr<TATensor<double,2>> I447, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I541)
   : ta0_(I447), ta1_(t2), ta2_(I541) { }
};

class Task315 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task315(std::shared_ptr<TATensor<double,4>> I541, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I541), ta1_(Gamma32), ta2_(t2) { }
};

class Task316 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task316(std::shared_ptr<TATensor<double,2>> I447, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I550)
   : ta0_(I447), ta1_(t2), ta2_(I550) { }
};

class Task317 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task317(std::shared_ptr<TATensor<double,4>> I550, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I550), ta1_(Gamma35), ta2_(t2) { }
};

class Task318 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task318(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I450)
   : ta0_(proj), ta1_(I450) { }
};

class Task319 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task319(std::shared_ptr<TATensor<double,2>> I450, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I451)
   : ta0_(I450), ta1_(t2), ta2_(I451) { }
};

class Task320 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task320(std::shared_ptr<TATensor<double,4>> I451, std::shared_ptr<TATensor<double,6>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I451), ta1_(Gamma2), ta2_(t2) { }
};

class Task321 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task321(std::shared_ptr<TATensor<double,2>> I450, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I556)
   : ta0_(I450), ta1_(t2), ta2_(I556) { }
};

class Task322 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task322(std::shared_ptr<TATensor<double,4>> I556, std::shared_ptr<TATensor<double,6>> Gamma37, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I556), ta1_(Gamma37), ta2_(t2) { }
};

class Task323 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task323(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I453)
   : ta0_(proj), ta1_(I453) { }
};

class Task324 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task324(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I454)
   : ta0_(I453), ta1_(t2), ta2_(I454) { }
};

class Task325 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task325(std::shared_ptr<TATensor<double,4>> I454, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I454), ta1_(Gamma3), ta2_(t2) { }
};

class Task326 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task326(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I565)
   : ta0_(I453), ta1_(t2), ta2_(I565) { }
};

class Task327 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task327(std::shared_ptr<TATensor<double,4>> I565, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I565), ta1_(Gamma35), ta2_(t2) { }
};

class Task328 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task328(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I568)
   : ta0_(I453), ta1_(t2), ta2_(I568) { }
};

class Task329 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task329(std::shared_ptr<TATensor<double,4>> I568, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I568), ta1_(Gamma32), ta2_(t2) { }
};

class Task330 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task330(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I607)
   : ta0_(I453), ta1_(t2), ta2_(I607) { }
};

class Task331 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task331(std::shared_ptr<TATensor<double,4>> I607, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I607), ta1_(Gamma35), ta2_(t2) { }
};

class Task332 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task332(std::shared_ptr<TATensor<double,2>> I453, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I610)
   : ta0_(I453), ta1_(t2), ta2_(I610) { }
};

class Task333 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task333(std::shared_ptr<TATensor<double,4>> I610, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I610), ta1_(Gamma35), ta2_(t2) { }
};

class Task334 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task334(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I456)
   : ta0_(proj), ta1_(I456) { }
};

class Task335 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task335(std::shared_ptr<TATensor<double,2>> I456, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I457)
   : ta0_(I456), ta1_(t2), ta2_(I457) { }
};

class Task336 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task336(std::shared_ptr<TATensor<double,4>> I457, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I457), ta1_(Gamma4), ta2_(t2) { }
};

class Task337 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task337(std::shared_ptr<TATensor<double,2>> I456, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I613)
   : ta0_(I456), ta1_(t2), ta2_(I613) { }
};

class Task338 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task338(std::shared_ptr<TATensor<double,4>> I613, std::shared_ptr<TATensor<double,6>> Gamma56, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I613), ta1_(Gamma56), ta2_(t2) { }
};

class Task339 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task339(std::shared_ptr<TATensor<double,2>> I456, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I616)
   : ta0_(I456), ta1_(t2), ta2_(I616) { }
};

class Task340 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task340(std::shared_ptr<TATensor<double,4>> I616, std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I616), ta1_(Gamma57), ta2_(t2) { }
};

class Task341 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task341(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I459)
   : ta0_(proj), ta1_(I459) { }
};

class Task342 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task342(std::shared_ptr<TATensor<double,2>> I459, std::shared_ptr<TATensor<double,8>> Gamma165, std::shared_ptr<TATensor<double,6>> I460)
   : ta0_(I459), ta1_(Gamma165), ta2_(I460) { }
};

class Task343 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task343(std::shared_ptr<TATensor<double,6>> I460, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I460), ta1_(t2) { }
};

class Task344 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task344(std::shared_ptr<TATensor<double,2>> I459, std::shared_ptr<TATensor<double,8>> Gamma218, std::shared_ptr<TATensor<double,6>> I619)
   : ta0_(I459), ta1_(Gamma218), ta2_(I619) { }
};

class Task345 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task345(std::shared_ptr<TATensor<double,6>> I619, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I619), ta1_(t2) { }
};

class Task346 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task346(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I462)
   : ta0_(proj), ta1_(I462) { }
};

class Task347 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task347(std::shared_ptr<TATensor<double,2>> I462, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I463)
   : ta0_(I462), ta1_(t2), ta2_(I463) { }
};

class Task348 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task348(std::shared_ptr<TATensor<double,4>> I463, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I463), ta1_(Gamma6), ta2_(t2) { }
};

class Task349 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task349(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I465)
   : ta0_(proj), ta1_(I465) { }
};


}
}
}
#endif
#endif

