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
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task700(std::shared_ptr<TATensor<std::complex<double>,6>> I1257, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1257), ta1_(t2), ta2_(v2) { }
};

class Task701 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task701(std::shared_ptr<TATensor<std::complex<double>,4>> I173, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma409, std::shared_ptr<TATensor<std::complex<double>,6>> I1260)
   : ta0_(I173), ta1_(Gamma409), ta2_(I1260) { }
};

class Task702 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task702(std::shared_ptr<TATensor<std::complex<double>,6>> I1260, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1260), ta1_(t2), ta2_(v2) { }
};

class Task703 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task703(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I194)
   : ta0_(proj), ta1_(I194) { }
};

class Task704 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task704(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> I195)
   : ta0_(I194), ta1_(Gamma2), ta2_(I195) { }
};

class Task705 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task705(std::shared_ptr<TATensor<std::complex<double>,4>> I195, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I195), ta1_(t2), ta2_(v2) { }
};

class Task706 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task706(std::shared_ptr<TATensor<std::complex<double>,4>> I195, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I195), ta1_(t2), ta2_(v2) { }
};

class Task707 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task707(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma412, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I194), ta1_(Gamma412), ta2_(t2) { }
};

class Task708 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task708(std::shared_ptr<TATensor<std::complex<double>,4>> I194, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma413, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I194), ta1_(Gamma413), ta2_(t2) { }
};

class Task709 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task709(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I773)
   : ta0_(proj), ta1_(I773) { }
};

class Task710 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task710(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I774)
   : ta0_(I773), ta1_(v2), ta2_(I774) { }
};

class Task711 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task711(std::shared_ptr<TATensor<std::complex<double>,4>> I774, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma2, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I774), ta1_(Gamma2), ta2_(t2) { }
};

class Task712 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task712(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task713 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task713(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task714 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task714(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task715 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task715(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I773), ta1_(t2), ta2_(v2) { }
};

class Task716 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task716(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I981)
   : ta0_(I773), ta1_(t2), ta2_(I981) { }
};

class Task717 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task717(std::shared_ptr<TATensor<std::complex<double>,4>> I981, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I981), ta1_(Gamma368), ta2_(v2) { }
};

class Task718 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,0>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task718(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,0>> Gamma424, std::shared_ptr<TATensor<std::complex<double>,4>> I1293)
   : ta0_(I773), ta1_(Gamma424), ta2_(I1293) { }
};

class Task719 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task719(std::shared_ptr<TATensor<std::complex<double>,4>> I1293, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1293), ta1_(t2) { }
};

class Task720 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,0>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task720(std::shared_ptr<TATensor<std::complex<double>,4>> I773, std::shared_ptr<TATensor<std::complex<double>,0>> Gamma426, std::shared_ptr<TATensor<std::complex<double>,4>> I1297)
   : ta0_(I773), ta1_(Gamma426), ta2_(I1297) { }
};

class Task721 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task721(std::shared_ptr<TATensor<std::complex<double>,4>> I1297, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1297), ta1_(t2) { }
};

class Task722 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task722(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1232)
   : ta0_(proj), ta1_(I1232) { }
};

class Task723 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task723(std::shared_ptr<TATensor<std::complex<double>,4>> I1232, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I1233)
   : ta0_(I1232), ta1_(t2), ta2_(I1233) { }
};

class Task724 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task724(std::shared_ptr<TATensor<std::complex<double>,4>> I1233, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1233), ta1_(Gamma368), ta2_(v2) { }
};

class Task725 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task725(std::shared_ptr<TATensor<std::complex<double>,4>> I1232, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma368, std::shared_ptr<TATensor<std::complex<double>,4>> I1266)
   : ta0_(I1232), ta1_(Gamma368), ta2_(I1266) { }
};

class Task726 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task726(std::shared_ptr<TATensor<std::complex<double>,4>> I1266, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1266), ta1_(t2), ta2_(v2) { }
};

class Task727 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task727(std::shared_ptr<TATensor<std::complex<double>,4>> I1232, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma432, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1232), ta1_(Gamma432), ta2_(t2) { }
};

class Task728 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task728(std::shared_ptr<TATensor<std::complex<double>,4>> I1232, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma433, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I1232), ta1_(Gamma433), ta2_(t2) { }
};

class Task729 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> tensor_;
    const bool reset_;
    void compute_() {
      if (reset_) tensor_->zero();
    }
  public:
    Task729(std::shared_ptr<TATensor<std::complex<double>,4>> t, const bool reset) : tensor_(t), reset_(reset) { }
};

class Task730 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task730(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1312)
   : ta0_(proj), ta1_(I1312) { }
};

class Task731 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task731(std::shared_ptr<TATensor<std::complex<double>,4>> I1312, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1312), ta1_(Gamma5), ta2_(h1) { }
};

class Task732 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task732(std::shared_ptr<TATensor<std::complex<double>,4>> I1312, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma81, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1312), ta1_(Gamma81), ta2_(v2) { }
};

class Task733 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task733(std::shared_ptr<TATensor<std::complex<double>,4>> I1312, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1312), ta1_(Gamma4), ta2_(v2) { }
};

class Task734 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task734(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1314)
   : ta0_(proj), ta1_(I1314) { }
};

class Task735 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task735(std::shared_ptr<TATensor<std::complex<double>,4>> I1314, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma27, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1314), ta1_(Gamma27), ta2_(h1) { }
};

class Task736 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task736(std::shared_ptr<TATensor<std::complex<double>,4>> I1314, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma24, std::shared_ptr<TATensor<std::complex<double>,4>> I1329)
   : ta0_(I1314), ta1_(Gamma24), ta2_(I1329) { }
};

class Task737 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task737(std::shared_ptr<TATensor<std::complex<double>,4>> I1329, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1329), ta1_(v2) { }
};

class Task738 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task738(std::shared_ptr<TATensor<std::complex<double>,4>> I1314, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma5, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1314), ta1_(Gamma5), ta2_(v2) { }
};

class Task739 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task739(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1316)
   : ta0_(proj), ta1_(I1316) { }
};

class Task740 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task740(std::shared_ptr<TATensor<std::complex<double>,4>> I1316, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I1316), ta1_(Gamma33), ta2_(h1) { }
};

class Task741 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task741(std::shared_ptr<TATensor<std::complex<double>,4>> I1316, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1316), ta1_(Gamma32), ta2_(v2) { }
};

class Task742 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task742(std::shared_ptr<TATensor<std::complex<double>,4>> I1316, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma31, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1316), ta1_(Gamma31), ta2_(v2) { }
};

class Task743 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task743(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1318)
   : ta0_(proj), ta1_(I1318) { }
};

class Task744 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task744(std::shared_ptr<TATensor<std::complex<double>,4>> I1318, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma0, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1318), ta1_(Gamma0), ta2_(v2) { }
};

class Task745 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task745(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1324)
   : ta0_(proj), ta1_(I1324) { }
};

class Task746 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task746(std::shared_ptr<TATensor<std::complex<double>,4>> I1324, std::shared_ptr<TATensor<std::complex<double>,2>> Gamma11, std::shared_ptr<TATensor<std::complex<double>,4>> I1325)
   : ta0_(I1324), ta1_(Gamma11), ta2_(I1325) { }
};

class Task747 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task747(std::shared_ptr<TATensor<std::complex<double>,4>> I1325, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1325), ta1_(v2) { }
};

class Task748 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task748(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I1340)
   : ta0_(proj), ta1_(I1340) { }
};

class Task749 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task749(std::shared_ptr<TATensor<std::complex<double>,4>> I1340, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I1340), ta1_(v2) { }
};


}
}
}
#endif
#endif

