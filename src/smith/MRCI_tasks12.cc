//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks12.cc
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

#include <src/smith/MRCI_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task550::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x5, a2, c1, x4") += (*ta1_)("x5, a2, c1, x4")
     + (*ta1_)("c1, a2, x5, x4") * (-2);
  madness::World::get_default().gop.fence();
}

void Task551::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x7, x6, x1, x0") * (*ta2_)("x7, a2, c1, x6");
  madness::World::get_default().gop.fence();
}

void Task552::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x7, a2, c1, x6") += (*ta1_)("x7, a2, c1, x6") * 0.5
     + (*ta1_)("c1, a2, x7, x6") * (-1);
  madness::World::get_default().gop.fence();
}

void Task553::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task554::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x3, x4, x2, x1") * (*ta2_)("x3, x5, a1, x4");
  madness::World::get_default().gop.fence();
}

void Task555::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta2_)("x3, c2");
  madness::World::get_default().gop.fence();
}

void Task556::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a2, c3, x4") * (*ta2_)("x3, c3, a2, a1");
  madness::World::get_default().gop.fence();
}

void Task557::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a1, c2, a3") * (*ta2_)("a3, x4, x3, c2") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task558::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a1, x5, x4");
  madness::World::get_default().gop.fence();
}

void Task559::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("x3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task560::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("x3, a1, a3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task561::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("c3, a2, x5, x4") * (*ta2_)("x3, c3, a2, a1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task562::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("x3, a1, a3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task563::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("x5, a3, c2, a1") * (*ta2_)("a3, x4, x3, c2") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task564::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("a1, x5, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task565::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a2, x4, x3") * (*ta2_)("a2, a1");
  madness::World::get_default().gop.fence();
}

void Task566::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a1, x4, a2") * (*ta2_)("a2, x3") * 2;
  madness::World::get_default().gop.fence();
}

void Task567::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a3, c2, a1") * (*ta2_)("a3, c2, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task568::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3") * (-0.5)
     + (*ta1_)("x4, x3, a3, c2") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task569::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a1, c2, a3") * (*ta2_)("a3, c2, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task570::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3")
     + (*ta1_)("x4, x3, a3, c2");
  madness::World::get_default().gop.fence();
}

void Task571::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a1, c3, a2") * (*ta2_)("x4, c3, a2, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task572::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a2, x4, a3") * (*ta2_)("a3, x3, a2, a1") * 2;
  madness::World::get_default().gop.fence();
}

void Task573::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1");
  madness::World::get_default().gop.fence();
}

void Task574::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task575::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task576::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("c2, a3, c4, a1") * (*ta2_)("x3, c4, a3, c2") * (-4);
  madness::World::get_default().gop.fence();
}

void Task577::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("c2, a1, c4, a3") * (*ta2_)("x3, c4, a3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task578::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("x3, a4, c2, a3") * (*ta2_)("a4, a1, a3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task579::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a1") += (*ta1_)("x3, a3, c2, a4") * (*ta2_)("a4, a1, a3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task580::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x4, a1, x3, c2") * (*ta2_)("c2, x3, x4, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task581::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x3, x4, x0, x2, x1") += (*ta1_)("x7, x6, x3, x5, x4, x0, x2, x1") * (*ta2_)("x7, x6, c2, x5") * (-1);
  madness::World::get_default().gop.fence();
}

void Task582::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x4, x5, x3, x0, x2, x1") * (*ta2_)("x4, x3, a1, x5");
  madness::World::get_default().gop.fence();
}

void Task583::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x4, x3, a1, x5") += (*ta1_)("c2, a1, c3, x5") * (*ta2_)("x4, c3, x3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task584::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, a1, c2, x6") * (*ta2_)("c2, x7, x0, x6, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task585::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x7, x0, x6, x2, x1") += (*ta1_)("x7, x0, x5, x6, x4, x3, x2, x1") * (*ta2_)("x5, c2, x4, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task586::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x7, x0, x6, x2, x1") += (*ta1_)("x7, x0, x3, x6, x5, x4, x2, x1") * (*ta2_)("x5, x4, x3, c2") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task587::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("c2, a1, x7, x6") * (*ta2_)("c2, x7, x6, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task588::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x7, x6, x0, x2, x1") += (*ta1_)("x7, x6, x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, c2, x4, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task589::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x7, x6, x0, x2, x1") += (*ta1_)("x7, x6, x5, x4, x3, x0, x2, x1") * (*ta2_)("x5, x4, x3, c2") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task590::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x0, x6, x5, x4, x3, x2, x1") * (*ta2_)("a1, x4, x3, x7, x6, x5");
  madness::World::get_default().gop.fence();
}

void Task591::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x4, x3, x7, x6, x5") += (*ta1_)("x7, a2, x6, x5") * (*ta2_)("a2, a1, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task592::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, a1, x4, x3") += (*ta1_)("a2, a1, x4, x3") * 0.5
     + (*ta1_)("x4, x3, a2, a1") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task593::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x4, x3, x7, x6, x5") += (*ta1_)("x7, a1, x6, a2") * (*ta2_)("a2, x5, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task594::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x4, x6, x5, x3, x0, x2, x1") * (*ta2_)("x4, x3, a1, x7, x6, x5");
  madness::World::get_default().gop.fence();
}

void Task595::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x4, x3, a1, x7, x6, x5") += (*ta1_)("x7, a2, x6, x5") * (*ta2_)("a2, x4, x3, a1") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task596::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x3, x6, x5, x4, x0, x2, x1") * (*ta2_)("x4, a1, x3, x7, x6, x5");
  madness::World::get_default().gop.fence();
}

void Task597::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x4, a1, x3, x7, x6, x5") += (*ta1_)("x7, a2, x6, x5") * (*ta2_)("x4, a1, a2, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task598::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x3, x4, x0, x2, x1") * (*ta2_)("x4, x3, x5, a1");
  madness::World::get_default().gop.fence();
}

void Task599::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x4, x3, x5, a1") += (*ta1_)("x5, a2, c3, a1") * (*ta2_)("x4, c3, a2, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

#endif
