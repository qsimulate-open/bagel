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
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task300(std::shared_ptr<TATensor<double,2>> I377, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I480)
   : ta0_(I377), ta1_(t2), ta2_(I480) { }
};

class Task301 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task301(std::shared_ptr<TATensor<double,4>> I480, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I480), ta1_(Gamma35), ta2_(t2) { }
};

class Task302 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task302(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I380)
   : ta0_(proj), ta1_(I380) { }
};

class Task303 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task303(std::shared_ptr<TATensor<double,2>> I380, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I381)
   : ta0_(I380), ta1_(t2), ta2_(I381) { }
};

class Task304 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task304(std::shared_ptr<TATensor<double,4>> I381, std::shared_ptr<TATensor<double,6>> Gamma2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I381), ta1_(Gamma2), ta2_(t2) { }
};

class Task305 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task305(std::shared_ptr<TATensor<double,2>> I380, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I486)
   : ta0_(I380), ta1_(t2), ta2_(I486) { }
};

class Task306 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task306(std::shared_ptr<TATensor<double,4>> I486, std::shared_ptr<TATensor<double,6>> Gamma37, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I486), ta1_(Gamma37), ta2_(t2) { }
};

class Task307 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task307(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I383)
   : ta0_(proj), ta1_(I383) { }
};

class Task308 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task308(std::shared_ptr<TATensor<double,2>> I383, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I384)
   : ta0_(I383), ta1_(t2), ta2_(I384) { }
};

class Task309 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task309(std::shared_ptr<TATensor<double,4>> I384, std::shared_ptr<TATensor<double,4>> Gamma3, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I384), ta1_(Gamma3), ta2_(t2) { }
};

class Task310 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task310(std::shared_ptr<TATensor<double,2>> I383, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I495)
   : ta0_(I383), ta1_(t2), ta2_(I495) { }
};

class Task311 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task311(std::shared_ptr<TATensor<double,4>> I495, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I495), ta1_(Gamma35), ta2_(t2) { }
};

class Task312 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task312(std::shared_ptr<TATensor<double,2>> I383, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I498)
   : ta0_(I383), ta1_(t2), ta2_(I498) { }
};

class Task313 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task313(std::shared_ptr<TATensor<double,4>> I498, std::shared_ptr<TATensor<double,4>> Gamma32, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I498), ta1_(Gamma32), ta2_(t2) { }
};

class Task314 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task314(std::shared_ptr<TATensor<double,2>> I383, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I537)
   : ta0_(I383), ta1_(t2), ta2_(I537) { }
};

class Task315 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task315(std::shared_ptr<TATensor<double,4>> I537, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I537), ta1_(Gamma35), ta2_(t2) { }
};

class Task316 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task316(std::shared_ptr<TATensor<double,2>> I383, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I540)
   : ta0_(I383), ta1_(t2), ta2_(I540) { }
};

class Task317 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task317(std::shared_ptr<TATensor<double,4>> I540, std::shared_ptr<TATensor<double,4>> Gamma35, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I540), ta1_(Gamma35), ta2_(t2) { }
};

class Task318 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task318(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I386)
   : ta0_(proj), ta1_(I386) { }
};

class Task319 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task319(std::shared_ptr<TATensor<double,2>> I386, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I387)
   : ta0_(I386), ta1_(t2), ta2_(I387) { }
};

class Task320 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task320(std::shared_ptr<TATensor<double,4>> I387, std::shared_ptr<TATensor<double,6>> Gamma4, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I387), ta1_(Gamma4), ta2_(t2) { }
};

class Task321 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task321(std::shared_ptr<TATensor<double,2>> I386, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I543)
   : ta0_(I386), ta1_(t2), ta2_(I543) { }
};

class Task322 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task322(std::shared_ptr<TATensor<double,4>> I543, std::shared_ptr<TATensor<double,6>> Gamma56, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I543), ta1_(Gamma56), ta2_(t2) { }
};

class Task323 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task323(std::shared_ptr<TATensor<double,2>> I386, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I546)
   : ta0_(I386), ta1_(t2), ta2_(I546) { }
};

class Task324 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task324(std::shared_ptr<TATensor<double,4>> I546, std::shared_ptr<TATensor<double,6>> Gamma57, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I546), ta1_(Gamma57), ta2_(t2) { }
};

class Task325 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task325(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I389)
   : ta0_(proj), ta1_(I389) { }
};

class Task326 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task326(std::shared_ptr<TATensor<double,2>> I389, std::shared_ptr<TATensor<double,8>> Gamma143, std::shared_ptr<TATensor<double,6>> I390)
   : ta0_(I389), ta1_(Gamma143), ta2_(I390) { }
};

class Task327 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task327(std::shared_ptr<TATensor<double,6>> I390, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I390), ta1_(t2) { }
};

class Task328 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task328(std::shared_ptr<TATensor<double,2>> I389, std::shared_ptr<TATensor<double,8>> Gamma196, std::shared_ptr<TATensor<double,6>> I549)
   : ta0_(I389), ta1_(Gamma196), ta2_(I549) { }
};

class Task329 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task329(std::shared_ptr<TATensor<double,6>> I549, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I549), ta1_(t2) { }
};

class Task330 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task330(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I392)
   : ta0_(proj), ta1_(I392) { }
};

class Task331 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task331(std::shared_ptr<TATensor<double,2>> I392, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I393)
   : ta0_(I392), ta1_(t2), ta2_(I393) { }
};

class Task332 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task332(std::shared_ptr<TATensor<double,4>> I393, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I393), ta1_(Gamma6), ta2_(t2) { }
};

class Task333 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task333(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I395)
   : ta0_(proj), ta1_(I395) { }
};

class Task334 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task334(std::shared_ptr<TATensor<double,2>> I395, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I396)
   : ta0_(I395), ta1_(t2), ta2_(I396) { }
};

class Task335 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task335(std::shared_ptr<TATensor<double,2>> I396, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I396), ta1_(Gamma7), ta2_(t2) { }
};

class Task336 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task336(std::shared_ptr<TATensor<double,2>> I395, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I399)
   : ta0_(I395), ta1_(t2), ta2_(I399) { }
};

class Task337 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task337(std::shared_ptr<TATensor<double,2>> I399, std::shared_ptr<TATensor<double,4>> Gamma7, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I399), ta1_(Gamma7), ta2_(t2) { }
};

class Task338 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task338(std::shared_ptr<TATensor<double,2>> I395, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I555)
   : ta0_(I395), ta1_(t2), ta2_(I555) { }
};

class Task339 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task339(std::shared_ptr<TATensor<double,2>> I555, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I555), ta1_(Gamma60), ta2_(t2) { }
};

class Task340 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task340(std::shared_ptr<TATensor<double,2>> I395, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> I558)
   : ta0_(I395), ta1_(t2), ta2_(I558) { }
};

class Task341 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task341(std::shared_ptr<TATensor<double,2>> I558, std::shared_ptr<TATensor<double,4>> Gamma60, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I558), ta1_(Gamma60), ta2_(t2) { }
};

class Task342 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task342(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I401)
   : ta0_(proj), ta1_(I401) { }
};

class Task343 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task343(std::shared_ptr<TATensor<double,2>> I401, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I402)
   : ta0_(I401), ta1_(t2), ta2_(I402) { }
};

class Task344 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task344(std::shared_ptr<TATensor<double,4>> I402, std::shared_ptr<TATensor<double,6>> Gamma9, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I402), ta1_(Gamma9), ta2_(t2) { }
};

class Task345 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task345(std::shared_ptr<TATensor<double,2>> I401, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I405)
   : ta0_(I401), ta1_(t2), ta2_(I405) { }
};

class Task346 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task346(std::shared_ptr<TATensor<double,4>> I405, std::shared_ptr<TATensor<double,6>> Gamma6, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I405), ta1_(Gamma6), ta2_(t2) { }
};

class Task347 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task347(std::shared_ptr<TATensor<double,2>> I401, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I561)
   : ta0_(I401), ta1_(t2), ta2_(I561) { }
};

class Task348 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task348(std::shared_ptr<TATensor<double,4>> I561, std::shared_ptr<TATensor<double,6>> Gamma59, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I561), ta1_(Gamma59), ta2_(t2) { }
};

class Task349 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    void compute_() override;
  public:
    Task349(std::shared_ptr<TATensor<double,2>> proj, std::shared_ptr<TATensor<double,2>> I407)
   : ta0_(proj), ta1_(I407) { }
};


}
}
}
#endif
#endif

