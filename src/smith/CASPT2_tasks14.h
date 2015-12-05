//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks14.h
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

#ifndef __SRC_SMITH_CASPT2_TASKS14_H
#define __SRC_SMITH_CASPT2_TASKS14_H

#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {
namespace CASPT2{

class Task650 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task650(std::shared_ptr<TATensor<double,4>> I881, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I882)
   : ta0_(I881), ta1_(t2), ta2_(I882) { }
};

class Task651 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task651(std::shared_ptr<TATensor<double,4>> I882, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I882), ta1_(t2), ta2_(f1) { }
};

class Task652 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task652(std::shared_ptr<TATensor<double,4>> I881, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I881), ta1_(v2), ta2_(t2) { }
};

class Task653 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task653(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma301, std::shared_ptr<TATensor<double,4>> I889)
   : ta0_(I768), ta1_(Gamma301), ta2_(I889) { }
};

class Task654 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task654(std::shared_ptr<TATensor<double,4>> I889, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I889), ta1_(t2) { }
};

class Task655 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task655(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma302, std::shared_ptr<TATensor<double,4>> I892)
   : ta0_(I768), ta1_(Gamma302), ta2_(I892) { }
};

class Task656 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task656(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I893)
   : ta0_(I892), ta1_(t2), ta2_(I893) { }
};

class Task657 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task657(std::shared_ptr<TATensor<double,4>> I893, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I893), ta1_(f1), ta2_(t2) { }
};

class Task658 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task658(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I897)
   : ta0_(I892), ta1_(t2), ta2_(I897) { }
};

class Task659 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task659(std::shared_ptr<TATensor<double,4>> I897, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I897), ta1_(f1), ta2_(t2) { }
};

class Task660 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task660(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I928)
   : ta0_(I892), ta1_(t2), ta2_(I928) { }
};

class Task661 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task661(std::shared_ptr<TATensor<double,4>> I928, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I928), ta1_(t2), ta2_(f1) { }
};

class Task662 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task662(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1055)
   : ta0_(I892), ta1_(t2), ta2_(I1055) { }
};

class Task663 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task663(std::shared_ptr<TATensor<double,4>> I1055, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I1055), ta1_(t2), e0_(e) { }
};

class Task664 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task664(std::shared_ptr<TATensor<double,4>> I1055, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1055), ta1_(f1), ta2_(t2) { }
};

class Task665 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task665(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I892), ta1_(v2), ta2_(t2) { }
};

class Task666 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task666(std::shared_ptr<TATensor<double,4>> I892, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I892), ta1_(v2), ta2_(t2) { }
};

class Task667 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task667(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma304, std::shared_ptr<TATensor<double,4>> I900)
   : ta0_(I768), ta1_(Gamma304), ta2_(I900) { }
};

class Task668 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task668(std::shared_ptr<TATensor<double,4>> I900, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I900), ta1_(t2) { }
};

class Task669 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task669(std::shared_ptr<TATensor<double,4>> I900, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I900), ta1_(t2) { }
};

class Task670 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    void compute_() override;
  public:
    Task670(std::shared_ptr<TATensor<double,4>> I900, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I900), ta1_(t2) { }
};

class Task671 : public Task {
  protected:
    std::shared_ptr<TATensor<double,1>> ta0_;
    std::shared_ptr<TATensor<double,5>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task671(std::shared_ptr<TATensor<double,1>> I768, std::shared_ptr<TATensor<double,5>> Gamma305, std::shared_ptr<TATensor<double,4>> I903)
   : ta0_(I768), ta1_(Gamma305), ta2_(I903) { }
};

class Task672 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task672(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I904)
   : ta0_(I903), ta1_(t2), ta2_(I904) { }
};

class Task673 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task673(std::shared_ptr<TATensor<double,4>> I904, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I904), ta1_(f1), ta2_(t2) { }
};

class Task674 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task674(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I908)
   : ta0_(I903), ta1_(t2), ta2_(I908) { }
};

class Task675 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task675(std::shared_ptr<TATensor<double,4>> I908, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I908), ta1_(f1), ta2_(t2) { }
};

class Task676 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task676(std::shared_ptr<TATensor<double,4>> I908, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I908), ta1_(f1), ta2_(t2) { }
};

class Task677 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task677(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I924)
   : ta0_(I903), ta1_(t2), ta2_(I924) { }
};

class Task678 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task678(std::shared_ptr<TATensor<double,4>> I924, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I924), ta1_(t2), ta2_(f1) { }
};

class Task679 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task679(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I947)
   : ta0_(I903), ta1_(t2), ta2_(I947) { }
};

class Task680 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task680(std::shared_ptr<TATensor<double,4>> I947, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I947), ta1_(f1), ta2_(t2) { }
};

class Task681 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task681(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I951)
   : ta0_(I903), ta1_(t2), ta2_(I951) { }
};

class Task682 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task682(std::shared_ptr<TATensor<double,4>> I951, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I951), ta1_(f1), ta2_(t2) { }
};

class Task683 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task683(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I958)
   : ta0_(I903), ta1_(t2), ta2_(I958) { }
};

class Task684 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task684(std::shared_ptr<TATensor<double,4>> I958, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I958), ta1_(f1), ta2_(t2) { }
};

class Task685 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task685(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I962)
   : ta0_(I903), ta1_(t2), ta2_(I962) { }
};

class Task686 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task686(std::shared_ptr<TATensor<double,4>> I962, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I962), ta1_(f1), ta2_(t2) { }
};

class Task687 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task687(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I978)
   : ta0_(I903), ta1_(t2), ta2_(I978) { }
};

class Task688 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task688(std::shared_ptr<TATensor<double,4>> I978, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I978), ta1_(t2), ta2_(f1) { }
};

class Task689 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,2>> ta2_;
    void compute_() override;
  public:
    Task689(std::shared_ptr<TATensor<double,4>> I978, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,2>> f1)
   : ta0_(I978), ta1_(t2), ta2_(f1) { }
};

class Task690 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task690(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1051)
   : ta0_(I903), ta1_(t2), ta2_(I1051) { }
};

class Task691 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task691(std::shared_ptr<TATensor<double,4>> I1051, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1051), ta1_(f1), ta2_(t2) { }
};

class Task692 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task692(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, std::shared_ptr<TATensor<double,4>> I1063)
   : ta0_(I903), ta1_(t2), ta2_(I1063) { }
};

class Task693 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task693(std::shared_ptr<TATensor<double,4>> I1063, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I1063), ta1_(t2), e0_(e) { }
};

class Task694 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,2>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task694(std::shared_ptr<TATensor<double,4>> I1063, std::shared_ptr<TATensor<double,2>> f1, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I1063), ta1_(f1), ta2_(t2) { }
};

class Task695 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task695(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I903), ta1_(t2), e0_(e) { }
};

class Task696 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    const double e0_;
    void compute_() override;
  public:
    Task696(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> t2, const double e)
   : ta0_(I903), ta1_(t2), e0_(e) { }
};

class Task697 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task697(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task698 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task698(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};

class Task699 : public Task {
  protected:
    std::shared_ptr<TATensor<double,4>> ta0_;
    std::shared_ptr<TATensor<double,4>> ta1_;
    std::shared_ptr<TATensor<double,4>> ta2_;
    void compute_() override;
  public:
    Task699(std::shared_ptr<TATensor<double,4>> I903, std::shared_ptr<TATensor<double,4>> v2, std::shared_ptr<TATensor<double,4>> t2)
   : ta0_(I903), ta1_(v2), ta2_(t2) { }
};


}
}
}
#endif
#endif

