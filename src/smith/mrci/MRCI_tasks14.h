//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks14.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#ifndef __SRC_SMITH_MRCI_TASKS14_H
#define __SRC_SMITH_MRCI_TASKS14_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task650 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task650(std::shared_ptr<TATensor<double,4>> I1135, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1135), ta1_(v2) { }
};

class Task651 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task651(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1137)
   : ta0_(I162), ta1_(t2), ta2_(I1137) { }
};

class Task652 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task652(std::shared_ptr<TATensor<double,4>> I1137, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1138)
   : ta0_(I1137), ta1_(Gamma12), ta2_(I1138) { }
};

class Task653 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task653(std::shared_ptr<TATensor<double,4>> I1138, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1138), ta1_(v2) { }
};

class Task654 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task654(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1146)
   : ta0_(I162), ta1_(t2), ta2_(I1146) { }
};

class Task655 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task655(std::shared_ptr<TATensor<double,4>> I1146, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1147)
   : ta0_(I1146), ta1_(Gamma12), ta2_(I1147) { }
};

class Task656 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task656(std::shared_ptr<TATensor<double,4>> I1147, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1147), ta1_(v2) { }
};

class Task657 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task657(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1149)
   : ta0_(I162), ta1_(t2), ta2_(I1149) { }
};

class Task658 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task658(std::shared_ptr<TATensor<double,4>> I1149, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> I1150)
   : ta0_(I1149), ta1_(Gamma12), ta2_(I1150) { }
};

class Task659 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task659(std::shared_ptr<TATensor<double,4>> I1150, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1150), ta1_(v2) { }
};

class Task660 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task660(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1158)
   : ta0_(I162), ta1_(v2), ta2_(I1158) { }
};

class Task661 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task661(std::shared_ptr<TATensor<double,4>> I1158, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1158), ta1_(Gamma12), ta2_(t2) { }
};

class Task662 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task662(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I1161)
   : ta0_(I162), ta1_(v2), ta2_(I1161) { }
};

class Task663 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task663(std::shared_ptr<TATensor<double,4>> I1161, std::shared_ptr<TATensor<double,2>> Gamma12, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1161), ta1_(Gamma12), ta2_(t2) { }
};

class Task664 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task664(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1164)
   : ta0_(I162), ta1_(t2), ta2_(I1164) { }
};

class Task665 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task665(std::shared_ptr<TATensor<double,4>> I1164, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1165)
   : ta0_(I1164), ta1_(Gamma29), ta2_(I1165) { }
};

class Task666 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task666(std::shared_ptr<TATensor<double,4>> I1165, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1165), ta1_(v2) { }
};

class Task667 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task667(std::shared_ptr<TATensor<double,4>> I1164, std::shared_ptr<TATensor<double,4>> Gamma18, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1164), ta1_(Gamma18), ta2_(v2) { }
};

class Task668 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task668(std::shared_ptr<TATensor<double,4>> I1164, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1164), ta1_(Gamma27), ta2_(v2) { }
};

class Task669 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task669(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1167)
   : ta0_(I162), ta1_(t2), ta2_(I1167) { }
};

class Task670 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task670(std::shared_ptr<TATensor<double,4>> I1167, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1168)
   : ta0_(I1167), ta1_(Gamma29), ta2_(I1168) { }
};

class Task671 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task671(std::shared_ptr<TATensor<double,4>> I1168, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1168), ta1_(v2) { }
};

class Task672 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task672(std::shared_ptr<TATensor<double,4>> I1167, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1167), ta1_(Gamma10), ta2_(v2) { }
};

class Task673 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task673(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1188)
   : ta0_(I162), ta1_(v2), ta2_(I1188) { }
};

class Task674 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task674(std::shared_ptr<TATensor<double,2>> I1188, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1189)
   : ta0_(I1188), ta1_(Gamma32), ta2_(I1189) { }
};

class Task675 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task675(std::shared_ptr<TATensor<double,4>> I1189, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1189), ta1_(t2) { }
};

class Task676 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task676(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1191)
   : ta0_(I162), ta1_(v2), ta2_(I1191) { }
};

class Task677 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task677(std::shared_ptr<TATensor<double,2>> I1191, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1192)
   : ta0_(I1191), ta1_(Gamma32), ta2_(I1192) { }
};

class Task678 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task678(std::shared_ptr<TATensor<double,4>> I1192, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1192), ta1_(t2) { }
};

class Task679 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task679(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1194)
   : ta0_(I162), ta1_(v2), ta2_(I1194) { }
};

class Task680 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task680(std::shared_ptr<TATensor<double,2>> I1194, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1195)
   : ta0_(I1194), ta1_(Gamma32), ta2_(I1195) { }
};

class Task681 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task681(std::shared_ptr<TATensor<double,4>> I1195, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1195), ta1_(t2) { }
};

class Task682 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task682(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1197)
   : ta0_(I162), ta1_(v2), ta2_(I1197) { }
};

class Task683 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task683(std::shared_ptr<TATensor<double,2>> I1197, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,4>> I1198)
   : ta0_(I1197), ta1_(Gamma32), ta2_(I1198) { }
};

class Task684 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task684(std::shared_ptr<TATensor<double,4>> I1198, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1198), ta1_(t2) { }
};

class Task685 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task685(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1200)
   : ta0_(I162), ta1_(t2), ta2_(I1200) { }
};

class Task686 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task686(std::shared_ptr<TATensor<double,4>> I1200, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1201)
   : ta0_(I1200), ta1_(Gamma29), ta2_(I1201) { }
};

class Task687 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task687(std::shared_ptr<TATensor<double,4>> I1201, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1201), ta1_(v2) { }
};

class Task688 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task688(std::shared_ptr<TATensor<double,4>> I1200, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1200), ta1_(Gamma10), ta2_(v2) { }
};

class Task689 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task689(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1203)
   : ta0_(I162), ta1_(t2), ta2_(I1203) { }
};

class Task690 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task690(std::shared_ptr<TATensor<double,4>> I1203, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I1204)
   : ta0_(I1203), ta1_(Gamma29), ta2_(I1204) { }
};

class Task691 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task691(std::shared_ptr<TATensor<double,4>> I1204, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1204), ta1_(v2) { }
};

class Task692 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task692(std::shared_ptr<TATensor<double,4>> I1203, std::shared_ptr<TATensor<double,4>> Gamma10, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I1203), ta1_(Gamma10), ta2_(v2) { }
};

class Task693 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task693(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1236)
   : ta0_(I162), ta1_(v2), ta2_(I1236) { }
};

class Task694 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task694(std::shared_ptr<TATensor<double,2>> I1236, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1236), ta1_(Gamma51), ta2_(t2) { }
};

class Task695 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task695(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,2>> I1239)
   : ta0_(I162), ta1_(v2), ta2_(I1239) { }
};

class Task696 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task696(std::shared_ptr<TATensor<double,2>> I1239, std::shared_ptr<TATensor<double,4>> Gamma51, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1239), ta1_(Gamma51), ta2_(t2) { }
};

class Task697 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task697(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task698 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task698(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};

class Task699 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task699(std::shared_ptr<TATensor<double,4>> I162, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I162), ta1_(t2), ta2_(v2) { }
};


}
}
}
#endif
#endif

