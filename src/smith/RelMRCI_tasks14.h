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
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task650(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1113)
   : ta0_(I134), ta1_(t2), ta2_(I1113) { }
};

class Task651 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task651(std::shared_ptr<TATensor<std::complex<double>,4>> I1113, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1114)
   : ta0_(I1113), ta1_(Gamma33), ta2_(I1114) { }
};

class Task652 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task652(std::shared_ptr<TATensor<std::complex<double>,4>> I1114, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1114), ta1_(v2) { }
};

class Task653 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task653(std::shared_ptr<TATensor<std::complex<double>,4>> I1113, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1113), ta1_(Gamma24), ta2_(v2) { }
};

class Task654 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task654(std::shared_ptr<TATensor<std::complex<double>,4>> I1113, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1113), ta1_(Gamma368), ta2_(v2) { }
};

class Task655 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task655(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1116)
   : ta0_(I134), ta1_(t2), ta2_(I1116) { }
};

class Task656 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task656(std::shared_ptr<TATensor<std::complex<double>,4>> I1116, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I1117)
   : ta0_(I1116), ta1_(Gamma33), ta2_(I1117) { }
};

class Task657 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task657(std::shared_ptr<TATensor<std::complex<double>,4>> I1117, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1117), ta1_(v2) { }
};

class Task658 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task658(std::shared_ptr<TATensor<std::complex<double>,4>> I1116, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma363, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1116), ta1_(Gamma363), ta2_(v2) { }
};

class Task659 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task659(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1203)
   : ta0_(I134), ta1_(t2), ta2_(I1203) { }
};

class Task660 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task660(std::shared_ptr<TATensor<std::complex<double>,4>> I1203, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1203), ta1_(Gamma32), ta2_(v2) { }
};

class Task661 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task661(std::shared_ptr<TATensor<std::complex<double>,4>> I1203, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma391, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1203), ta1_(Gamma391), ta2_(v2) { }
};

class Task662 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task662(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1209)
   : ta0_(I134), ta1_(v2), ta2_(I1209) { }
};

class Task663 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task663(std::shared_ptr<TATensor<std::complex<double>,4>> I1209, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1209), ta1_(Gamma33), ta2_(t2) { }
};

class Task664 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task664(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1212)
   : ta0_(I134), ta1_(v2), ta2_(I1212) { }
};

class Task665 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task665(std::shared_ptr<TATensor<std::complex<double>,4>> I1212, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1212), ta1_(Gamma33), ta2_(t2) { }
};

class Task666 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task666(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1215)
   : ta0_(I134), ta1_(v2), ta2_(I1215) { }
};

class Task667 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task667(std::shared_ptr<TATensor<std::complex<double>,4>> I1215, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1215), ta1_(Gamma368), ta2_(t2) { }
};

class Task668 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task668(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1218)
   : ta0_(I134), ta1_(v2), ta2_(I1218) { }
};

class Task669 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task669(std::shared_ptr<TATensor<std::complex<double>,4>> I1218, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1218), ta1_(Gamma33), ta2_(t2) { }
};

class Task670 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task670(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma428, std::shared_ptr<TATensor<std::complex<double>,4>> I1301)
   : ta0_(I134), ta1_(Gamma428), ta2_(I1301) { }
};

class Task671 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task671(std::shared_ptr<TATensor<std::complex<double>,4>> I1301, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1301), ta1_(t2) { }
};

class Task672 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task672(std::shared_ptr<TATensor<std::complex<double>,4>> I134, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma430, std::shared_ptr<TATensor<std::complex<double>,4>> I1305)
   : ta0_(I134), ta1_(Gamma430), ta2_(I1305) { }
};

class Task673 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task673(std::shared_ptr<TATensor<std::complex<double>,4>> I1305, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1305), ta1_(t2) { }
};

class Task674 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task674(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I173)
   : ta0_(proj), ta1_(I173) { }
};

class Task675 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task675(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,2>> h1, std::shared_ptr<TATensor<std::complex<double>,4>> I174)
   : ta0_(I173), ta1_(h1), ta2_(I174) { }
};

class Task676 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task676(std::shared_ptr<TATensor<std::complex<double>,4>> I174, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I174), ta1_(Gamma32), ta2_(t2) { }
};

class Task677 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task677(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> I177)
   : ta0_(I173), ta1_(Gamma33), ta2_(I177) { }
};

class Task678 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task678(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I177), ta1_(t2), ta2_(h1) { }
};

class Task679 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task679(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I177), ta1_(t2), ta2_(h1) { }
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
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task681(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I177), ta1_(t2), ta2_(v2) { }
};

class Task682 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task682(std::shared_ptr<TATensor<std::complex<double>,4>> I177, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I177), ta1_(t2), ta2_(v2) { }
};

class Task683 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task683(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I1221)
   : ta0_(I173), ta1_(t2), ta2_(I1221) { }
};

class Task684 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task684(std::shared_ptr<TATensor<std::complex<double>,6>> I1221, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma396, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1221), ta1_(Gamma396), ta2_(v2) { }
};

class Task685 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task685(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I1224)
   : ta0_(I173), ta1_(t2), ta2_(I1224) { }
};

class Task686 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task686(std::shared_ptr<TATensor<std::complex<double>,6>> I1224, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma397, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1224), ta1_(Gamma397), ta2_(v2) { }
};

class Task687 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task687(std::shared_ptr<TATensor<std::complex<double>,6>> I1224, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma236, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1224), ta1_(Gamma236), ta2_(v2) { }
};

class Task688 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task688(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I1230)
   : ta0_(I173), ta1_(v2), ta2_(I1230) { }
};

class Task689 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task689(std::shared_ptr<TATensor<std::complex<double>,4>> I1230, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma336, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1230), ta1_(Gamma336), ta2_(t2) { }
};

class Task690 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task690(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1236)
   : ta0_(I173), ta1_(t2), ta2_(I1236) { }
};

class Task691 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task691(std::shared_ptr<TATensor<std::complex<double>,4>> I1236, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma391, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1236), ta1_(Gamma391), ta2_(v2) { }
};

class Task692 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task692(std::shared_ptr<TATensor<std::complex<double>,4>> I1236, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1236), ta1_(Gamma32), ta2_(v2) { }
};

class Task693 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task693(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> I1242)
   : ta0_(I173), ta1_(Gamma368), ta2_(I1242) { }
};

class Task694 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task694(std::shared_ptr<TATensor<std::complex<double>,4>> I1242, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1242), ta1_(t2), ta2_(v2) { }
};

class Task695 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task695(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma391, std::shared_ptr<TATensor<std::complex<double>,6>> I1254)
   : ta0_(I173), ta1_(Gamma391), ta2_(I1254) { }
};

class Task696 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task696(std::shared_ptr<TATensor<std::complex<double>,6>> I1254, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1255)
   : ta0_(I1254), ta1_(t2), ta2_(I1255) { }
};

class Task697 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task697(std::shared_ptr<TATensor<std::complex<double>,4>> I1255, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1255), ta1_(v2) { }
};

class Task698 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task698(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,6>> I1257)
   : ta0_(I173), ta1_(Gamma32), ta2_(I1257) { }
};

class Task699 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task699(std::shared_ptr<TATensor<std::complex<double>,6>> I1257, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1257), ta1_(t2), ta2_(v2) { }
};


}
}
}
#endif
#endif

