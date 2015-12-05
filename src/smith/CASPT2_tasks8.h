//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks8.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS8_H
#define __SRC_SMITH_CASPT2_TASKS8_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task350 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task350(std::shared_ptr<TATensor<double,2>> I465, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I466)
   : ta0_(I465), ta1_(t2), ta2_(I466) { }
};

class Task351 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task351(std::shared_ptr<TATensor<double,2>> I466, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I466), ta1_(Gamma7), ta2_(t2) { }
};

class Task352 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task352(std::shared_ptr<TATensor<double,2>> I465, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I469)
   : ta0_(I465), ta1_(t2), ta2_(I469) { }
};

class Task353 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task353(std::shared_ptr<TATensor<double,2>> I469, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I469), ta1_(Gamma7), ta2_(t2) { }
};

class Task354 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task354(std::shared_ptr<TATensor<double,2>> I465, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I625)
   : ta0_(I465), ta1_(t2), ta2_(I625) { }
};

class Task355 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task355(std::shared_ptr<TATensor<double,2>> I625, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I625), ta1_(Gamma60), ta2_(t2) { }
};

class Task356 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task356(std::shared_ptr<TATensor<double,2>> I465, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I628)
   : ta0_(I465), ta1_(t2), ta2_(I628) { }
};

class Task357 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task357(std::shared_ptr<TATensor<double,2>> I628, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I628), ta1_(Gamma60), ta2_(t2) { }
};

class Task358 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task358(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I471)
   : ta0_(proj), ta1_(I471) { }
};

class Task359 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task359(std::shared_ptr<TATensor<double,2>> I471, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I472)
   : ta0_(I471), ta1_(t2), ta2_(I472) { }
};

class Task360 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task360(std::shared_ptr<TATensor<double,4>> I472, std::shared_ptr<TATensor<double,6>> Gamma9, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I472), ta1_(Gamma9), ta2_(t2) { }
};

class Task361 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task361(std::shared_ptr<TATensor<double,2>> I471, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I475)
   : ta0_(I471), ta1_(t2), ta2_(I475) { }
};

class Task362 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task362(std::shared_ptr<TATensor<double,4>> I475, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I475), ta1_(Gamma6), ta2_(t2) { }
};

class Task363 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task363(std::shared_ptr<TATensor<double,2>> I471, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I631)
   : ta0_(I471), ta1_(t2), ta2_(I631) { }
};

class Task364 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task364(std::shared_ptr<TATensor<double,4>> I631, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I631), ta1_(Gamma59), ta2_(t2) { }
};

class Task365 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task365(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I477)
   : ta0_(proj), ta1_(I477) { }
};

class Task366 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task366(std::shared_ptr<TATensor<double,2>> I477, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I478)
   : ta0_(I477), ta1_(t2), ta2_(I478) { }
};

class Task367 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task367(std::shared_ptr<TATensor<double,4>> I478, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I478), ta1_(Gamma3), ta2_(t2) { }
};

class Task368 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task368(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I480)
   : ta0_(proj), ta1_(I480) { }
};

class Task369 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task369(std::shared_ptr<TATensor<double,2>> I480, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I481)
   : ta0_(I480), ta1_(t2), ta2_(I481) { }
};

class Task370 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task370(std::shared_ptr<TATensor<double,2>> I481, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I481), ta1_(Gamma12), ta2_(t2) { }
};

class Task371 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task371(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I483)
   : ta0_(proj), ta1_(I483) { }
};

class Task372 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task372(std::shared_ptr<TATensor<double,2>> I483, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I484)
   : ta0_(I483), ta1_(t2), ta2_(I484) { }
};

class Task373 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task373(std::shared_ptr<TATensor<double,2>> I484, std::shared_ptr<TATensor<double,4>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I484), ta1_(Gamma12), ta2_(t2) { }
};

class Task374 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task374(std::shared_ptr<TATensor<double,2>> I483, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I640)
   : ta0_(I483), ta1_(t2), ta2_(I640) { }
};

class Task375 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task375(std::shared_ptr<TATensor<double,2>> I640, std::shared_ptr<TATensor<double,2>> Gamma38, std::shared_ptr<TATensor<double,4>> I641)
   : ta0_(I640), ta1_(Gamma38), ta2_(I641) { }
};

class Task376 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task376(std::shared_ptr<TATensor<double,4>> I641, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I641), ta1_(t2) { }
};

class Task377 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task377(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I486)
   : ta0_(proj), ta1_(I486) { }
};

class Task378 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task378(std::shared_ptr<TATensor<double,2>> I486, std::shared_ptr<TATensor<double,4>> Gamma174, std::shared_ptr<TATensor<double,2>> I487)
   : ta0_(I486), ta1_(Gamma174), ta2_(I487) { }
};

class Task379 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task379(std::shared_ptr<TATensor<double,2>> I487, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I487), ta1_(t2) { }
};

class Task380 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task380(std::shared_ptr<TATensor<double,2>> I487, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I487), ta1_(t2) { }
};

class Task381 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task381(std::shared_ptr<TATensor<double,2>> I486, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,2>> I696)
   : ta0_(I486), ta1_(Gamma60), ta2_(I696) { }
};

class Task382 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task382(std::shared_ptr<TATensor<double,2>> I696, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I696), ta1_(t2) { }
};

class Task383 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task383(std::shared_ptr<TATensor<double,2>> I696, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I696), ta1_(t2) { }
};

class Task384 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task384(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I492)
   : ta0_(proj), ta1_(I492) { }
};

class Task385 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task385(std::shared_ptr<TATensor<double,2>> I492, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I493)
   : ta0_(I492), ta1_(t2), ta2_(I493) { }
};

class Task386 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task386(std::shared_ptr<TATensor<double,4>> I493, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I493), ta1_(Gamma16), ta2_(t2) { }
};

class Task387 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task387(std::shared_ptr<TATensor<double,2>> I492, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I496)
   : ta0_(I492), ta1_(t2), ta2_(I496) { }
};

class Task388 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task388(std::shared_ptr<TATensor<double,4>> I496, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I496), ta1_(Gamma16), ta2_(t2) { }
};

class Task389 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task389(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I498)
   : ta0_(proj), ta1_(I498) { }
};

class Task390 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task390(std::shared_ptr<TATensor<double,2>> I498, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I499)
   : ta0_(I498), ta1_(t2), ta2_(I499) { }
};

class Task391 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task391(std::shared_ptr<TATensor<double,4>> I499, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I499), ta1_(Gamma16), ta2_(t2) { }
};

class Task392 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task392(std::shared_ptr<TATensor<double,2>> I498, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I505)
   : ta0_(I498), ta1_(t2), ta2_(I505) { }
};

class Task393 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task393(std::shared_ptr<TATensor<double,4>> I505, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I505), ta1_(Gamma16), ta2_(t2) { }
};

class Task394 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task394(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I501)
   : ta0_(proj), ta1_(I501) { }
};

class Task395 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task395(std::shared_ptr<TATensor<double,2>> I501, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I502)
   : ta0_(I501), ta1_(t2), ta2_(I502) { }
};

class Task396 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task396(std::shared_ptr<TATensor<double,4>> I502, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I502), ta1_(Gamma16), ta2_(t2) { }
};

class Task397 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task397(std::shared_ptr<TATensor<double,2>> I501, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I508)
   : ta0_(I501), ta1_(t2), ta2_(I508) { }
};

class Task398 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task398(std::shared_ptr<TATensor<double,4>> I508, std::shared_ptr<TATensor<double,2>> Gamma16, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I508), ta1_(Gamma16), ta2_(t2) { }
};

class Task399 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task399(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I510)
   : ta0_(proj), ta1_(I510) { }
};


}
}
}
#endif
#endif

