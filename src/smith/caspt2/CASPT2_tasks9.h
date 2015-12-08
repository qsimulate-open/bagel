//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks9.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS9_H
#define __SRC_SMITH_CASPT2_TASKS9_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task400 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task400(std::shared_ptr<TATensor<double,2>> I510, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I511)
   : ta0_(I510), ta1_(t2), ta2_(I511) { }
};

class Task401 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task401(std::shared_ptr<TATensor<double,4>> I511, std::shared_ptr<TATensor<double,4>> Gamma22, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I511), ta1_(Gamma22), ta2_(t2) { }
};

class Task402 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task402(std::shared_ptr<TATensor<double,4>> I511, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I511), ta1_(Gamma12), ta2_(t2) { }
};

class Task403 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task403(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I513)
   : ta0_(proj), ta1_(I513) { }
};

class Task404 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task404(std::shared_ptr<TATensor<double,2>> I513, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I514)
   : ta0_(I513), ta1_(t2), ta2_(I514) { }
};

class Task405 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task405(std::shared_ptr<TATensor<double,4>> I514, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> I515)
   : ta0_(I514), ta1_(Gamma12), ta2_(I515) { }
};

class Task406 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task406(std::shared_ptr<TATensor<double,4>> I515, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I515), ta1_(t2) { }
};

class Task407 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task407(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I522)
   : ta0_(proj), ta1_(I522) { }
};

class Task408 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task408(std::shared_ptr<TATensor<double,2>> I522, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,2>> I523)
   : ta0_(I522), ta1_(Gamma16), ta2_(I523) { }
};

class Task409 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task409(std::shared_ptr<TATensor<double,2>> I523, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I523), ta1_(t2) { }
};

class Task410 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task410(std::shared_ptr<TATensor<double,2>> I523, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I523), ta1_(t2) { }
};

class Task411 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task411(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I528)
   : ta0_(proj), ta1_(I528) { }
};

class Task412 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task412(std::shared_ptr<TATensor<double,2>> I528, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I529)
   : ta0_(I528), ta1_(t2), ta2_(I529) { }
};

class Task413 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task413(std::shared_ptr<TATensor<double,4>> I529, std::shared_ptr<TATensor<double,6>> Gamma28, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I529), ta1_(Gamma28), ta2_(t2) { }
};

class Task414 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task414(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I531)
   : ta0_(proj), ta1_(I531) { }
};

class Task415 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task415(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I532)
   : ta0_(I531), ta1_(t2), ta2_(I532) { }
};

class Task416 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task416(std::shared_ptr<TATensor<double,4>> I532, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I532), ta1_(Gamma29), ta2_(t2) { }
};

class Task417 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task417(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I535)
   : ta0_(I531), ta1_(t2), ta2_(I535) { }
};

class Task418 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task418(std::shared_ptr<TATensor<double,4>> I535, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I535), ta1_(Gamma7), ta2_(t2) { }
};

class Task419 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task419(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I574)
   : ta0_(I531), ta1_(t2), ta2_(I574) { }
};

class Task420 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task420(std::shared_ptr<TATensor<double,4>> I574, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I574), ta1_(Gamma7), ta2_(t2) { }
};

class Task421 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task421(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I577)
   : ta0_(I531), ta1_(t2), ta2_(I577) { }
};

class Task422 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task422(std::shared_ptr<TATensor<double,4>> I577, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I577), ta1_(Gamma7), ta2_(t2) { }
};

class Task423 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task423(std::shared_ptr<TATensor<double,2>> I531, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I726)
   : ta0_(I531), ta1_(t2), ta2_(I726) { }
};

class Task424 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task424(std::shared_ptr<TATensor<double,4>> I726, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I726), ta1_(Gamma60), ta2_(t2) { }
};

class Task425 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task425(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I543)
   : ta0_(proj), ta1_(I543) { }
};

class Task426 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task426(std::shared_ptr<TATensor<double,2>> I543, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I544)
   : ta0_(I543), ta1_(t2), ta2_(I544) { }
};

class Task427 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task427(std::shared_ptr<TATensor<double,4>> I544, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I544), ta1_(Gamma32), ta2_(t2) { }
};

class Task428 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task428(std::shared_ptr<TATensor<double,2>> I543, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I553)
   : ta0_(I543), ta1_(t2), ta2_(I553) { }
};

class Task429 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task429(std::shared_ptr<TATensor<double,4>> I553, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I553), ta1_(Gamma35), ta2_(t2) { }
};

class Task430 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task430(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I558)
   : ta0_(proj), ta1_(I558) { }
};

class Task431 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task431(std::shared_ptr<TATensor<double,2>> I558, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I559)
   : ta0_(I558), ta1_(t2), ta2_(I559) { }
};

class Task432 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task432(std::shared_ptr<TATensor<double,2>> I559, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I559), ta1_(Gamma38), ta2_(t2) { }
};

class Task433 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task433(std::shared_ptr<TATensor<double,2>> I558, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I562)
   : ta0_(I558), ta1_(t2), ta2_(I562) { }
};

class Task434 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task434(std::shared_ptr<TATensor<double,2>> I562, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I562), ta1_(Gamma38), ta2_(t2) { }
};

class Task435 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task435(std::shared_ptr<TATensor<double,2>> I558, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I601)
   : ta0_(I558), ta1_(t2), ta2_(I601) { }
};

class Task436 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task436(std::shared_ptr<TATensor<double,2>> I601, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I601), ta1_(Gamma38), ta2_(t2) { }
};

class Task437 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task437(std::shared_ptr<TATensor<double,2>> I558, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I604)
   : ta0_(I558), ta1_(t2), ta2_(I604) { }
};

class Task438 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task438(std::shared_ptr<TATensor<double,2>> I604, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I604), ta1_(Gamma38), ta2_(t2) { }
};

class Task439 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task439(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I570)
   : ta0_(proj), ta1_(I570) { }
};

class Task440 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task440(std::shared_ptr<TATensor<double,2>> I570, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I571)
   : ta0_(I570), ta1_(t2), ta2_(I571) { }
};

class Task441 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task441(std::shared_ptr<TATensor<double,4>> I571, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I571), ta1_(Gamma6), ta2_(t2) { }
};

class Task442 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task442(std::shared_ptr<TATensor<double,2>> I570, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I723)
   : ta0_(I570), ta1_(t2), ta2_(I723) { }
};

class Task443 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task443(std::shared_ptr<TATensor<double,4>> I723, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I723), ta1_(Gamma59), ta2_(t2) { }
};

class Task444 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task444(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I582)
   : ta0_(proj), ta1_(I582) { }
};

class Task445 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task445(std::shared_ptr<TATensor<double,2>> I582, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I583)
   : ta0_(I582), ta1_(t2), ta2_(I583) { }
};

class Task446 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task446(std::shared_ptr<TATensor<double,4>> I583, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I583), ta1_(Gamma35), ta2_(t2) { }
};

class Task447 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task447(std::shared_ptr<TATensor<double,2>> I582, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I592)
   : ta0_(I582), ta1_(t2), ta2_(I592) { }
};

class Task448 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task448(std::shared_ptr<TATensor<double,4>> I592, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I592), ta1_(Gamma35), ta2_(t2) { }
};

class Task449 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task449(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I585)
   : ta0_(proj), ta1_(I585) { }
};


}
}
}
#endif
#endif

