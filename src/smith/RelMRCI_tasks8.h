//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks8.h
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

#ifndef __SRC_SMITH_RelMRCI_TASKS8_H
#define __SRC_SMITH_RelMRCI_TASKS8_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class Task350 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task350(std::shared_ptr<TATensor<std::complex<double>,4>> I558, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma183, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I558), ta1_(Gamma183), ta2_(v2) { }
};

class Task351 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task351(std::shared_ptr<TATensor<std::complex<double>,4>> I558, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma81, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I558), ta1_(Gamma81), ta2_(v2) { }
};

class Task352 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task352(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I561)
   : ta0_(I63), ta1_(t2), ta2_(I561) { }
};

class Task353 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task353(std::shared_ptr<TATensor<std::complex<double>,4>> I561, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma183, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I561), ta1_(Gamma183), ta2_(v2) { }
};

class Task354 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task354(std::shared_ptr<TATensor<std::complex<double>,4>> I561, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma81, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I561), ta1_(Gamma81), ta2_(v2) { }
};

class Task355 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task355(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I588)
   : ta0_(I63), ta1_(t2), ta2_(I588) { }
};

class Task356 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task356(std::shared_ptr<TATensor<std::complex<double>,6>> I588, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma193, std::shared_ptr<TATensor<std::complex<double>,4>> I589)
   : ta0_(I588), ta1_(Gamma193), ta2_(I589) { }
};

class Task357 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task357(std::shared_ptr<TATensor<std::complex<double>,4>> I589, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I589), ta1_(v2) { }
};

class Task358 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task358(std::shared_ptr<TATensor<std::complex<double>,6>> I588, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I588), ta1_(Gamma4), ta2_(v2) { }
};

class Task359 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task359(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma193, std::shared_ptr<TATensor<std::complex<double>,6>> I591)
   : ta0_(I63), ta1_(Gamma193), ta2_(I591) { }
};

class Task360 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task360(std::shared_ptr<TATensor<std::complex<double>,6>> I591, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I592)
   : ta0_(I591), ta1_(t2), ta2_(I592) { }
};

class Task361 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task361(std::shared_ptr<TATensor<std::complex<double>,4>> I592, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I592), ta1_(v2) { }
};

class Task362 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task362(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma4, std::shared_ptr<TATensor<std::complex<double>,6>> I597)
   : ta0_(I63), ta1_(Gamma4), ta2_(I597) { }
};

class Task363 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task363(std::shared_ptr<TATensor<std::complex<double>,6>> I597, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I597), ta1_(t2), ta2_(v2) { }
};

class Task364 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task364(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,6>> I618)
   : ta0_(I63), ta1_(t2), ta2_(I618) { }
};

class Task365 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task365(std::shared_ptr<TATensor<std::complex<double>,6>> I618, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma203, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I618), ta1_(Gamma203), ta2_(v2) { }
};

class Task366 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,8>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task366(std::shared_ptr<TATensor<std::complex<double>,6>> I618, std::shared_ptr<TATensor<std::complex<double>,8>> Gamma204, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I618), ta1_(Gamma204), ta2_(v2) { }
};

class Task367 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task367(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I624)
   : ta0_(I63), ta1_(v2), ta2_(I624) { }
};

class Task368 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task368(std::shared_ptr<TATensor<std::complex<double>,4>> I624, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I624), ta1_(Gamma26), ta2_(t2) { }
};

class Task369 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task369(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> v2, std::shared_ptr<TATensor<std::complex<double>,4>> I627)
   : ta0_(I63), ta1_(v2), ta2_(I627) { }
};

class Task370 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task370(std::shared_ptr<TATensor<std::complex<double>,4>> I627, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I627), ta1_(Gamma26), ta2_(t2) { }
};

class Task371 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task371(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I666)
   : ta0_(I63), ta1_(t2), ta2_(I666) { }
};

class Task372 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task372(std::shared_ptr<TATensor<std::complex<double>,4>> I666, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma193, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I666), ta1_(Gamma193), ta2_(v2) { }
};

class Task373 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task373(std::shared_ptr<TATensor<std::complex<double>,4>> I666, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I666), ta1_(Gamma26), ta2_(v2) { }
};

class Task374 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task374(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I669)
   : ta0_(I63), ta1_(t2), ta2_(I669) { }
};

class Task375 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task375(std::shared_ptr<TATensor<std::complex<double>,4>> I669, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma193, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I669), ta1_(Gamma193), ta2_(v2) { }
};

class Task376 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task376(std::shared_ptr<TATensor<std::complex<double>,4>> I669, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma26, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I669), ta1_(Gamma26), ta2_(v2) { }
};

class Task377 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta2_;
    void compute_() override;
  public:
    Task377(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma229, std::shared_ptr<TATensor<std::complex<double>,6>> I696)
   : ta0_(I63), ta1_(Gamma229), ta2_(I696) { }
};

class Task378 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,6>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task378(std::shared_ptr<TATensor<std::complex<double>,6>> I696, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I696), ta1_(t2), ta2_(v2) { }
};

class Task379 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task379(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma420, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I63), ta1_(Gamma420), ta2_(t2) { }
};

class Task380 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task380(std::shared_ptr<TATensor<std::complex<double>,4>> I63, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma421, std::shared_ptr<TATensor<std::complex<double>,4>> t2)
   : ta0_(I63), ta1_(Gamma421), ta2_(t2) { }
};

class Task381 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task381(std::shared_ptr<TATensor<std::complex<double>,4>> proj, std::shared_ptr<TATensor<std::complex<double>,4>> I93)
   : ta0_(proj), ta1_(I93) { }
};

class Task382 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task382(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma31, std::shared_ptr<TATensor<std::complex<double>,4>> I94)
   : ta0_(I93), ta1_(Gamma31), ta2_(I94) { }
};

class Task383 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task383(std::shared_ptr<TATensor<std::complex<double>,4>> I94, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I94), ta1_(t2), ta2_(h1) { }
};

class Task384 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task384(std::shared_ptr<TATensor<std::complex<double>,4>> I94, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I94), ta1_(t2), ta2_(v2) { }
};

class Task385 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task385(std::shared_ptr<TATensor<std::complex<double>,4>> I94, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I94), ta1_(t2), ta2_(v2) { }
};

class Task386 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task386(std::shared_ptr<TATensor<std::complex<double>,4>> I94, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I94), ta1_(t2), ta2_(v2) { }
};

class Task387 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,6>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task387(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,6>> Gamma32, std::shared_ptr<TATensor<std::complex<double>,4>> I97)
   : ta0_(I93), ta1_(Gamma32), ta2_(I97) { }
};

class Task388 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task388(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I97), ta1_(t2), ta2_(h1) { }
};

class Task389 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task389(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I97), ta1_(t2), ta2_(h1) { }
};

class Task390 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task390(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I736)
   : ta0_(I97), ta1_(t2), ta2_(I736) { }
};

class Task391 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task391(std::shared_ptr<TATensor<std::complex<double>,4>> I736, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I736), ta1_(v2) { }
};

class Task392 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task392(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> I739)
   : ta0_(I97), ta1_(t2), ta2_(I739) { }
};

class Task393 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    void compute_() override;
  public:
    Task393(std::shared_ptr<TATensor<std::complex<double>,4>> I739, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I739), ta1_(v2) { }
};

class Task394 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task394(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task395 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task395(std::shared_ptr<TATensor<std::complex<double>,4>> I97, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task396 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,4>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task396(std::shared_ptr<TATensor<std::complex<double>,4>> I93, std::shared_ptr<TATensor<std::complex<double>,4>> Gamma33, std::shared_ptr<TATensor<std::complex<double>,2>> I100)
   : ta0_(I93), ta1_(Gamma33), ta2_(I100) { }
};

class Task397 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task397(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I100), ta1_(t2), ta2_(h1) { }
};

class Task398 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,2>> ta2_;
    void compute_() override;
  public:
    Task398(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,2>> h1)
   : ta0_(I100), ta1_(t2), ta2_(h1) { }
};

class Task399 : public Task {
  protected:
    std::shared_ptr<TATensor<std::complex<double>,2>> ta0_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta1_;
    std::shared_ptr<TATensor<std::complex<double>,4>> ta2_;
    void compute_() override;
  public:
    Task399(std::shared_ptr<TATensor<std::complex<double>,2>> I100, std::shared_ptr<TATensor<std::complex<double>,4>> t2, std::shared_ptr<TATensor<std::complex<double>,4>> v2)
   : ta0_(I100), ta1_(t2), ta2_(v2) { }
};


}
}
}
#endif
#endif

