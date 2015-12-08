//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks15.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS15_H
#define __SRC_SMITH_RelMRCI_TASKS15_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task700 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task700(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> I195)
   : ta0_(I194), ta1_(Gamma2), ta2_(I195) { }
};

class Task701 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task701(std::shared_ptr<TATensor<std::complex<double>,4>> I195, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I195), ta1_(t2), ta2_(v2) { }
};

class Task702 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task702(std::shared_ptr<TATensor<std::complex<double>,4>> I195, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I195), ta1_(t2), ta2_(v2) { }
};

class Task703 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task703(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma409, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I194), ta1_(Gamma409), ta2_(t2) { }
};

class Task704 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task704(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma410, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I194), ta1_(Gamma410), ta2_(t2) { }
};

class Task705 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task705(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I773)
   : ta0_(proj), ta1_(I773) { }
};

class Task706 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task706(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I774)
   : ta0_(I773), ta1_(v2), ta2_(I774) { }
};

class Task707 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task707(std::shared_ptr<TATensor<std::complex<double>,4>> I774, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I774), ta1_(Gamma2), ta2_(t2) { }
};

class Task708 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task708(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task709 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task709(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task710 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task710(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I977)
   : ta0_(I773), ta1_(t2), ta2_(I977) { }
};

class Task711 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task711(std::shared_ptr<TATensor<std::complex<double>,4>> I977, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I977), ta1_(Gamma368), ta2_(v2) { }
};

class Task712 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,0>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task712(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,0>> Gamma421, std::shared_ptr<TATensor<std::complex<double>,4>> I1280)
   : ta0_(I773), ta1_(Gamma421), ta2_(I1280) { }
};

class Task713 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task713(std::shared_ptr<TATensor<std::complex<double>,4>> I1280, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1280), ta1_(t2) { }
};

class Task714 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,0>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task714(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,0>> Gamma423, std::shared_ptr<TATensor<std::complex<double>,4>> I1284)
   : ta0_(I773), ta1_(Gamma423), ta2_(I1284) { }
};

class Task715 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task715(std::shared_ptr<TATensor<std::complex<double>,4>> I1284, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1284), ta1_(t2) { }
};

class Task716 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task716(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1222)
   : ta0_(proj), ta1_(I1222) { }
};

class Task717 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task717(std::shared_ptr<TATensor<std::complex<double>,4>> I1222, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1223)
   : ta0_(I1222), ta1_(t2), ta2_(I1223) { }
};

class Task718 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task718(std::shared_ptr<TATensor<std::complex<double>,4>> I1223, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1223), ta1_(Gamma368), ta2_(v2) { }
};

class Task719 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task719(std::shared_ptr<TATensor<std::complex<double>,4>> I1222, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma429, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1222), ta1_(Gamma429), ta2_(t2) { }
};

class Task720 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task720(std::shared_ptr<TATensor<std::complex<double>,4>> I1222, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma430, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1222), ta1_(Gamma430), ta2_(t2) { }
};

class Task721 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task721(std::shared_ptr<TATensor<std::complex<double>,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task722 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task722(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1299)
   : ta0_(proj), ta1_(I1299) { }
};

class Task723 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task723(std::shared_ptr<TATensor<std::complex<double>,4>> I1299, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1299), ta1_(Gamma5), ta2_(h1) { }
};

class Task724 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task724(std::shared_ptr<TATensor<std::complex<double>,4>> I1299, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma81, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1299), ta1_(Gamma81), ta2_(v2) { }
};

class Task725 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task725(std::shared_ptr<TATensor<std::complex<double>,4>> I1299, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1299), ta1_(Gamma4), ta2_(v2) { }
};

class Task726 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task726(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1301)
   : ta0_(proj), ta1_(I1301) { }
};

class Task727 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task727(std::shared_ptr<TATensor<std::complex<double>,4>> I1301, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1301), ta1_(Gamma27), ta2_(h1) { }
};

class Task728 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task728(std::shared_ptr<TATensor<std::complex<double>,4>> I1301, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> I1316)
   : ta0_(I1301), ta1_(Gamma24), ta2_(I1316) { }
};

class Task729 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task729(std::shared_ptr<TATensor<std::complex<double>,4>> I1316, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1316), ta1_(v2) { }
};

class Task730 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task730(std::shared_ptr<TATensor<std::complex<double>,4>> I1301, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1301), ta1_(Gamma5), ta2_(v2) { }
};

class Task731 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task731(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1303)
   : ta0_(proj), ta1_(I1303) { }
};

class Task732 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task732(std::shared_ptr<TATensor<std::complex<double>,4>> I1303, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1303), ta1_(Gamma33), ta2_(h1) { }
};

class Task733 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task733(std::shared_ptr<TATensor<std::complex<double>,4>> I1303, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1303), ta1_(Gamma32), ta2_(v2) { }
};

class Task734 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task734(std::shared_ptr<TATensor<std::complex<double>,4>> I1303, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma31, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1303), ta1_(Gamma31), ta2_(v2) { }
};

class Task735 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task735(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1305)
   : ta0_(proj), ta1_(I1305) { }
};

class Task736 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task736(std::shared_ptr<TATensor<std::complex<double>,4>> I1305, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1305), ta1_(Gamma0), ta2_(v2) { }
};

class Task737 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task737(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1311)
   : ta0_(proj), ta1_(I1311) { }
};

class Task738 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task738(std::shared_ptr<TATensor<std::complex<double>,4>> I1311, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I1312)
   : ta0_(I1311), ta1_(Gamma11), ta2_(I1312) { }
};

class Task739 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task739(std::shared_ptr<TATensor<std::complex<double>,4>> I1312, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1312), ta1_(v2) { }
};

class Task740 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task740(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1327)
   : ta0_(proj), ta1_(I1327) { }
};

class Task741 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task741(std::shared_ptr<TATensor<std::complex<double>,4>> I1327, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1327), ta1_(v2) { }
};

class Task742 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task742(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1329)
   : ta0_(proj), ta1_(I1329) { }
};

class Task743 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task743(std::shared_ptr<TATensor<std::complex<double>,4>> I1329, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,4>> I1330)
   : ta0_(I1329), ta1_(Gamma27), ta2_(I1330) { }
};

class Task744 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task744(std::shared_ptr<TATensor<std::complex<double>,4>> I1330, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1330), ta1_(v2) { }
};

class Task745 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task745(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1333)
   : ta0_(proj), ta1_(I1333) { }
};

class Task746 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task746(std::shared_ptr<TATensor<std::complex<double>,4>> I1333, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1333), ta1_(Gamma33), ta2_(v2) { }
};

class Task747 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task747(std::shared_ptr<TATensor<std::complex<double>,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task748 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task748(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1335)
   : ta0_(proj), ta1_(I1335) { }
};

class Task749 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task749(std::shared_ptr<TATensor<std::complex<double>,4>> I1335, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1335), ta1_(Gamma0), ta2_(t2) { }
};


}
}
}
#endif
#endif

