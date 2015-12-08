//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks14.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS14_H
#define __SRC_SMITH_RelMRCI_TASKS14_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task650 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task650(std::shared_ptr<TATensor<std::complex<double>,4>> I1110, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1110), ta1_(v2) { }
};

class Task651 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task651(std::shared_ptr<TATensor<std::complex<double>,4>> I1109, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1109), ta1_(Gamma24), ta2_(v2) { }
};

class Task652 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task652(std::shared_ptr<TATensor<std::complex<double>,4>> I1109, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1109), ta1_(Gamma368), ta2_(v2) { }
};

class Task653 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task653(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1112)
   : ta0_(I134), ta1_(t2), ta2_(I1112) { }
};

class Task654 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task654(std::shared_ptr<TATensor<std::complex<double>,4>> I1112, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1113)
   : ta0_(I1112), ta1_(Gamma33), ta2_(I1113) { }
};

class Task655 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task655(std::shared_ptr<TATensor<std::complex<double>,4>> I1113, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1113), ta1_(v2) { }
};

class Task656 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task656(std::shared_ptr<TATensor<std::complex<double>,4>> I1112, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma363, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1112), ta1_(Gamma363), ta2_(v2) { }
};

class Task657 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task657(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1193)
   : ta0_(I134), ta1_(t2), ta2_(I1193) { }
};

class Task658 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task658(std::shared_ptr<TATensor<std::complex<double>,4>> I1193, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1193), ta1_(Gamma32), ta2_(v2) { }
};

class Task659 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task659(std::shared_ptr<TATensor<std::complex<double>,4>> I1193, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma389, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1193), ta1_(Gamma389), ta2_(v2) { }
};

class Task660 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task660(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1199)
   : ta0_(I134), ta1_(v2), ta2_(I1199) { }
};

class Task661 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task661(std::shared_ptr<TATensor<std::complex<double>,4>> I1199, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1199), ta1_(Gamma33), ta2_(t2) { }
};

class Task662 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task662(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1202)
   : ta0_(I134), ta1_(v2), ta2_(I1202) { }
};

class Task663 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task663(std::shared_ptr<TATensor<std::complex<double>,4>> I1202, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1202), ta1_(Gamma33), ta2_(t2) { }
};

class Task664 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task664(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1205)
   : ta0_(I134), ta1_(v2), ta2_(I1205) { }
};

class Task665 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task665(std::shared_ptr<TATensor<std::complex<double>,4>> I1205, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1205), ta1_(Gamma368), ta2_(t2) { }
};

class Task666 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task666(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1208)
   : ta0_(I134), ta1_(v2), ta2_(I1208) { }
};

class Task667 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task667(std::shared_ptr<TATensor<std::complex<double>,4>> I1208, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1208), ta1_(Gamma33), ta2_(t2) { }
};

class Task668 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task668(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma425, std::shared_ptr<TATensor<std::complex<double>,4>> I1288)
   : ta0_(I134), ta1_(Gamma425), ta2_(I1288) { }
};

class Task669 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task669(std::shared_ptr<TATensor<std::complex<double>,4>> I1288, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1288), ta1_(t2) { }
};

class Task670 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task670(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma427, std::shared_ptr<TATensor<std::complex<double>,4>> I1292)
   : ta0_(I134), ta1_(Gamma427), ta2_(I1292) { }
};

class Task671 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task671(std::shared_ptr<TATensor<std::complex<double>,4>> I1292, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1292), ta1_(t2) { }
};

class Task672 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task672(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I173)
   : ta0_(proj), ta1_(I173) { }
};

class Task673 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task673(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I174)
   : ta0_(I173), ta1_(h1), ta2_(I174) { }
};

class Task674 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task674(std::shared_ptr<TATensor<std::complex<double>,4>> I174, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I174), ta1_(Gamma32), ta2_(t2) { }
};

class Task675 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task675(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I177)
   : ta0_(I173), ta1_(Gamma33), ta2_(I177) { }
};

class Task676 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task676(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I177), ta1_(t2), ta2_(h1) { }
};

class Task677 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task677(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I177), ta1_(t2), ta2_(h1) { }
};

class Task678 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task678(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I177), ta1_(t2), ta2_(v2) { }
};

class Task679 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task679(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I177), ta1_(t2), ta2_(v2) { }
};

class Task680 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task680(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I177), ta1_(t2), ta2_(v2) { }
};

class Task681 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task681(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I1211)
   : ta0_(I173), ta1_(t2), ta2_(I1211) { }
};

class Task682 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task682(std::shared_ptr<TATensor<std::complex<double>,6>> I1211, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma394, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1211), ta1_(Gamma394), ta2_(v2) { }
};

class Task683 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task683(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I1214)
   : ta0_(I173), ta1_(t2), ta2_(I1214) { }
};

class Task684 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task684(std::shared_ptr<TATensor<std::complex<double>,6>> I1214, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma395, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1214), ta1_(Gamma395), ta2_(v2) { }
};

class Task685 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task685(std::shared_ptr<TATensor<std::complex<double>,6>> I1214, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma236, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1214), ta1_(Gamma236), ta2_(v2) { }
};

class Task686 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task686(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1220)
   : ta0_(I173), ta1_(v2), ta2_(I1220) { }
};

class Task687 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task687(std::shared_ptr<TATensor<std::complex<double>,4>> I1220, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma336, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1220), ta1_(Gamma336), ta2_(t2) { }
};

class Task688 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task688(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1226)
   : ta0_(I173), ta1_(t2), ta2_(I1226) { }
};

class Task689 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task689(std::shared_ptr<TATensor<std::complex<double>,4>> I1226, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma389, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1226), ta1_(Gamma389), ta2_(v2) { }
};

class Task690 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task690(std::shared_ptr<TATensor<std::complex<double>,4>> I1226, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1226), ta1_(Gamma32), ta2_(v2) { }
};

class Task691 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task691(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> I1232)
   : ta0_(I173), ta1_(Gamma368), ta2_(I1232) { }
};

class Task692 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task692(std::shared_ptr<TATensor<std::complex<double>,4>> I1232, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1232), ta1_(t2), ta2_(v2) { }
};

class Task693 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task693(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma389, std::shared_ptr<TATensor<std::complex<double>,6>> I1244)
   : ta0_(I173), ta1_(Gamma389), ta2_(I1244) { }
};

class Task694 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task694(std::shared_ptr<TATensor<std::complex<double>,6>> I1244, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1245)
   : ta0_(I1244), ta1_(t2), ta2_(I1245) { }
};

class Task695 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task695(std::shared_ptr<TATensor<std::complex<double>,4>> I1245, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1245), ta1_(v2) { }
};

class Task696 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task696(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,6>> I1247)
   : ta0_(I173), ta1_(Gamma32), ta2_(I1247) { }
};

class Task697 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task697(std::shared_ptr<TATensor<std::complex<double>,6>> I1247, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1247), ta1_(t2), ta2_(v2) { }
};

class Task698 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task698(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma407, std::shared_ptr<TATensor<std::complex<double>,6>> I1250)
   : ta0_(I173), ta1_(Gamma407), ta2_(I1250) { }
};

class Task699 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task699(std::shared_ptr<TATensor<std::complex<double>,6>> I1250, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1250), ta1_(t2), ta2_(v2) { }
};


}
}
}
#endif
#endif

