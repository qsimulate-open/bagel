//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks8.h
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

#ifndef __SRC_SMITH_MRCI_TASKS8_H
#define __SRC_SMITH_MRCI_TASKS8_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace MRCI{

class Task350 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task350(std::shared_ptr<TATensor<double,4>> I76, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I76), ta1_(t2), ta2_(v2) { }
};

class Task351 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task351(std::shared_ptr<TATensor<double,4>> I76, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I76), ta1_(t2), ta2_(v2) { }
};

class Task352 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task352(std::shared_ptr<TATensor<double,4>> I76, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I76), ta1_(t2), ta2_(v2) { }
};

class Task353 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task353(std::shared_ptr<TATensor<double,4>> I76, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I76), ta1_(t2), ta2_(v2) { }
};

class Task354 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task354(std::shared_ptr<TATensor<double,4>> I76, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I76), ta1_(t2), ta2_(v2) { }
};

class Task355 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task355(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma5, std::shared_ptr<TATensor<double,4>> I79)
   : ta0_(I72), ta1_(Gamma5), ta2_(I79) { }
};

class Task356 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task356(std::shared_ptr<TATensor<double,4>> I79, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I79), ta1_(t2), ta2_(h1) { }
};

class Task357 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task357(std::shared_ptr<TATensor<double,4>> I79, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I79), ta1_(t2), ta2_(v2) { }
};

class Task358 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task358(std::shared_ptr<TATensor<double,4>> I79, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I79), ta1_(t2), ta2_(v2) { }
};

class Task359 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task359(std::shared_ptr<TATensor<double,4>> I79, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I79), ta1_(t2), ta2_(v2) { }
};

class Task360 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task360(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma27, std::shared_ptr<TATensor<double,4>> I82)
   : ta0_(I72), ta1_(Gamma27), ta2_(I82) { }
};

class Task361 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task361(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I82), ta1_(t2), ta2_(h1) { }
};

class Task362 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task362(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I82), ta1_(t2), ta2_(h1) { }
};

class Task363 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task363(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I82), ta1_(t2), ta2_(h1) { }
};

class Task364 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task364(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task365 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task365(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task366 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task366(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I823)
   : ta0_(I82), ta1_(t2), ta2_(I823) { }
};

class Task367 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task367(std::shared_ptr<TATensor<double,4>> I823, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I823), ta1_(v2) { }
};

class Task368 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task368(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task369 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task369(std::shared_ptr<TATensor<double,4>> I82, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I82), ta1_(t2), ta2_(v2) { }
};

class Task370 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task370(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> Gamma29, std::shared_ptr<TATensor<double,4>> I88)
   : ta0_(I72), ta1_(Gamma29), ta2_(I88) { }
};

class Task371 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task371(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I88), ta1_(t2), ta2_(h1) { }
};

class Task372 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task372(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I88), ta1_(t2), ta2_(h1) { }
};

class Task373 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task373(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I88), ta1_(t2), ta2_(h1) { }
};

class Task374 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task374(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task375 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task375(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task376 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task376(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task377 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task377(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I772)
   : ta0_(I88), ta1_(t2), ta2_(I772) { }
};

class Task378 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task378(std::shared_ptr<TATensor<double,4>> I772, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I772), ta1_(v2) { }
};

class Task379 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task379(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I775)
   : ta0_(I88), ta1_(t2), ta2_(I775) { }
};

class Task380 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task380(std::shared_ptr<TATensor<double,4>> I775, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I775), ta1_(v2) { }
};

class Task381 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task381(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task382 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task382(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task383 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task383(std::shared_ptr<TATensor<double,4>> I88, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I88), ta1_(t2), ta2_(v2) { }
};

class Task384 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task384(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,2>> h1, std::shared_ptr<TATensor<double,4>> I94)
   : ta0_(I72), ta1_(h1), ta2_(I94) { }
};

class Task385 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task385(std::shared_ptr<TATensor<double,4>> I94, std::shared_ptr<TATensor<double,6>> Gamma31, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I94), ta1_(Gamma31), ta2_(t2) { }
};

class Task386 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task386(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,2>> Gamma32, std::shared_ptr<TATensor<double,2>> I97)
   : ta0_(I72), ta1_(Gamma32), ta2_(I97) { }
};

class Task387 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task387(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I97), ta1_(t2), ta2_(h1) { }
};

class Task388 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task388(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> h1)
   : ta0_(I97), ta1_(t2), ta2_(h1) { }
};

class Task389 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task389(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task390 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task390(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task391 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task391(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task392 : public Task {
  protected:
    std::shared_ptr<TATensor<double,2>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task392(std::shared_ptr<TATensor<double,2>> I97, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> v2)
   : ta0_(I97), ta1_(t2), ta2_(v2) { }
};

class Task393 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task393(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I654)
   : ta0_(I72), ta1_(v2), ta2_(I654) { }
};

class Task394 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,6>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task394(std::shared_ptr<TATensor<double,6>> I654, std::shared_ptr<TATensor<double,6>> Gamma215, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I654), ta1_(Gamma215), ta2_(t2) { }
};

class Task395 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task395(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I657)
   : ta0_(I72), ta1_(v2), ta2_(I657) { }
};

class Task396 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task396(std::shared_ptr<TATensor<double,6>> I657, std::shared_ptr<TATensor<double,8>> Gamma216, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I657), ta1_(Gamma216), ta2_(t2) { }
};

class Task397 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,6>> ta2_;
    void compute_() override;
  public:
    Task397(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,6>> I660)
   : ta0_(I72), ta1_(v2), ta2_(I660) { }
};

class Task398 : public Task {
  protected:
    std::shared_ptr<TATensor<double,6>> ta0_;
    std::shared_ptr<TATensor<double,8>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task398(std::shared_ptr<TATensor<double,6>> I660, std::shared_ptr<TATensor<double,8>> Gamma217, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I660), ta1_(Gamma217), ta2_(t2) { }
};

class Task399 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task399(std::shared_ptr<TATensor<double,4>> I72, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> I663)
   : ta0_(I72), ta1_(v2), ta2_(I663) { }
};


}
}
}
#endif
#endif

