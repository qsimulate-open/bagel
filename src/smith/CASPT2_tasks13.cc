//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks13.cc
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

#include <src/smith/CASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task600::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("a2, x3") * (-2);
  madness::World::get_default().gop.fence();
}

void Task601::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x2, x0, x1") * (*ta2_)("x0, x3, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task602::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("x3, x2, c3, x1") * (*ta2_)("x0, c3");
  madness::World::get_default().gop.fence();
}

void Task603::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3") += (*ta1_)("c1, a2") * (*ta2_)("c1, a2, c3, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task604::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("x3, x2, c1, x1") * (*ta2_)("x0, c1");
  madness::World::get_default().gop.fence();
}

void Task605::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c1") += (*ta1_)("c3, a2") * (*ta2_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task606::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("x3, a2, c1, x2") * (*ta2_)("x0, a2, c1, x1");
  madness::World::get_default().gop.fence();
}

void Task607::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a2, c1, x1") += (*ta1_)("c3, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task608::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("c3, a2, x3, x2") * (*ta2_)("x0, c3, a2, x1");
  madness::World::get_default().gop.fence();
}

void Task609::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, a2, x1") += (*ta1_)("c1, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task610::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("c1, a2, x3, x2") * (*ta2_)("x0, a2, c1, x1");
  madness::World::get_default().gop.fence();
}

void Task611::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a2, c1, x1") += (*ta1_)("c3, x1") * (*ta2_)("c1, a2, c3, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task612::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("c1, x0, x1, a2") * (*ta2_)("c1, a2, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task613::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3, x2, x1") += (*ta1_)("c1, x0") * (*ta2_)("x3, x2, c1, x1") * 2;
  madness::World::get_default().gop.fence();
}

void Task614::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("ci0") += (*ta1_)("ci0, x0, x3") * (*ta2_)("x0, x3");
  madness::World::get_default().gop.fence();
}

void Task615::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3") += (*ta1_)("c3, a2, c1, x3") * (*ta1_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task616::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x3") += (*ta1_)("c1, a2, c3, x3") * (*ta1_)("c1, a2, c3, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task617::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("ci0") += (*ta1_)("ci0, x0, x1") * (*ta2_)("x0, x1");
  madness::World::get_default().gop.fence();
}

void Task618::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c4, a2, c3, x1") * (*ta2_)("x0, c3, a2, c4");
  madness::World::get_default().gop.fence();
}

void Task619::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, a2, c4") += (*ta1_)("c1, c4") * (*ta2_)("c1, a2, c3, x0") * (-4);
  madness::World::get_default().gop.fence();
}

void Task620::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c4, x1") * (*ta2_)("x0, c3, a2, c4");
  madness::World::get_default().gop.fence();
}

void Task621::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, a2, c4") += (*ta1_)("c1, c4") * (*ta2_)("c1, a2, c3, x0") * 2;
  madness::World::get_default().gop.fence();
}

void Task622::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c4, a2, c1, x1") * (*ta2_)("x0, a2, c1, c4");
  madness::World::get_default().gop.fence();
}

void Task623::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a2, c1, c4") += (*ta1_)("c3, c4") * (*ta2_)("c1, a2, c3, x0") * 2;
  madness::World::get_default().gop.fence();
}

void Task624::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c3, a4, c1, x1") * (*ta2_)("x0, c3, c1, a4");
  madness::World::get_default().gop.fence();
}

void Task625::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, c1, a4") += (*ta1_)("a4, a2") * (*ta2_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task626::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c4, x1") * (*ta2_)("x0, a2, c1, c4");
  madness::World::get_default().gop.fence();
}

void Task627::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a2, c1, c4") += (*ta1_)("c3, c4") * (*ta2_)("c1, a2, c3, x0") * (-4);
  madness::World::get_default().gop.fence();
}

void Task628::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("x0, c3, c1, a4");
  madness::World::get_default().gop.fence();
}

void Task629::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, c1, a4") += (*ta1_)("a4, a2") * (*ta2_)("c1, a2, c3, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task630::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("a4, x1") * (*ta2_)("x0, a4");
  madness::World::get_default().gop.fence();
}

void Task631::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a4") += (*ta1_)("c1, a4, c3, a2") * (*ta1_)("c1, a2, c3, x0") * (-4);
  madness::World::get_default().gop.fence();
}

void Task632::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, a4") += (*ta1_)("c1, a2, c3, a4") * (*ta1_)("c1, a2, c3, x0") * 8;
  madness::World::get_default().gop.fence();
}

void Task633::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("x0, a2") * (*ta2_)("a2, x1");
  madness::World::get_default().gop.fence();
}

void Task634::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta1_)("c1, a2, c3, a4") * (-4);
  madness::World::get_default().gop.fence();
}

void Task635::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("x0, a4") * (*ta2_)("a4, x1");
  madness::World::get_default().gop.fence();
}

void Task636::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, a4") * 8;
  madness::World::get_default().gop.fence();
}

void Task637::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c1, x1") * (*ta1_)("c1, a2, c3, x0") * e0_ * 2;
  madness::World::get_default().gop.fence();
}

void Task638::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, x0") * e0_ * (-4);
  madness::World::get_default().gop.fence();
}

void Task639::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c3, x1, c1, a2") * (*ta2_)("c1, a2, c3, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task640::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c1, x1, c3, a2") * (*ta2_)("c1, a2, c3, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task641::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c3, x0, c1, a2") * (*ta2_)("c1, a2, c3, x1") * 4;
  madness::World::get_default().gop.fence();
}

void Task642::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("c1, a2, c3, x1") * (-2);
  madness::World::get_default().gop.fence();
}

void Task643::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x1, x0, x2") * (*ta2_)("x0, x1, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task644::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1, x3, x2") += (*ta1_)("x3, a2, c3, x2") * (*ta2_)("x0, c3, a2, x1");
  madness::World::get_default().gop.fence();
}

void Task645::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c3, a2, x1") += (*ta1_)("c1, x1") * (*ta2_)("c1, a2, c3, x0") * 2;
  madness::World::get_default().gop.fence();
}

void Task646::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, x1, x3, x2") += (*ta1_)("c2, x0, x1, a1") * (*ta2_)("x3, a1, c2, x2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task647::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x1, x3, x2, x0") * (*ta2_)("x1, x0, x2, x5, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task648::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, x0, x2, x5, x4, x3") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("x1, c2, x0, x2");
  madness::World::get_default().gop.fence();
}

void Task649::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, c2, x0, x2") += (*ta1_)("x2, a1") * (*ta2_)("x0, a1, c2, x1") * (-2);
  madness::World::get_default().gop.fence();
}

#endif
