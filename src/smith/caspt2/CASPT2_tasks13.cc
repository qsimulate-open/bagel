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

#include <src/smith/caspt2/CASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task600::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x0, x1") * (*ta2_)("x0, x1");
}

void Task601::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c4, a2, c3, x1") * (*ta2_)("x0, c3, a2, c4");
}

void Task602::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c3, a2, c4") += (*ta1_)("c1, c4") * (*ta2_)("c1, a2, c3, x0") * (-4);
}

void Task603::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c4, x1") * (*ta2_)("x0, c3, a2, c4");
}

void Task604::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c3, a2, c4") += (*ta1_)("c1, c4") * (*ta2_)("c1, a2, c3, x0") * 2;
}

void Task605::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c4, a2, c1, x1") * (*ta2_)("x0, a2, c1, c4");
}

void Task606::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a2, c1, c4") += (*ta1_)("c3, c4") * (*ta2_)("c1, a2, c3, x0") * 2;
}

void Task607::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, a4, c1, x1") * (*ta2_)("x0, c3, c1, a4");
}

void Task608::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c3, c1, a4") += (*ta1_)("a4, a2") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task609::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c4, x1") * (*ta2_)("x0, a2, c1, c4");
}

void Task610::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a2, c1, c4") += (*ta1_)("c3, c4") * (*ta2_)("c1, a2, c3, x0") * (-4);
}

void Task611::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("x0, c3, c1, a4");
}

void Task612::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c3, c1, a4") += (*ta1_)("a4, a2") * (*ta2_)("c1, a2, c3, x0") * 4;
}

void Task613::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("a4, x1") * (*ta2_)("x0, a4");
}

void Task614::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a4, c3, a2") * (*ta1_)("c1, a2, c3, x0") * (-4);
}

void Task615::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a2, c3, a4") * (*ta1_)("c1, a2, c3, x0") * 8;
}

void Task616::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x0, a2") * (*ta2_)("a2, x1");
}

void Task617::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta1_)("c1, a2, c3, a4") * (-4);
}

void Task618::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x0, a4") * (*ta2_)("a4, x1");
}

void Task619::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, a4") * 8;
}

void Task620::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c1, x1") * (*ta1_)("c1, a2, c3, x0") * e0_ * 2;
}

void Task621::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, x0") * e0_ * (-4);
}

void Task622::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, x1, c1, a2") * (*ta2_)("c1, a2, c3, x0") * 4;
}

void Task623::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, x1, c3, a2") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task624::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, x0, c1, a2") * (*ta2_)("c1, a2, c3, x1") * 4;
}

void Task625::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("c1, a2, c3, x1") * (-2);
}

void Task626::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x1, x0, x2") * (*ta2_)("x0, x1, x3, x2");
}

void Task627::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1, x3, x2") += (*ta1_)("x3, a2, c3, x2") * (*ta2_)("x0, c3, a2, x1");
}

void Task628::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c3, a2, x1") += (*ta1_)("c1, x1") * (*ta2_)("c1, a2, c3, x0") * 2;
}

void Task629::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1, x3, x2") += (*ta1_)("c2, x0, x1, a1") * (*ta2_)("x3, a1, c2, x2") * (-1);
}

void Task630::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x1, x3, x2, x0") * (*ta2_)("x1, x0, x2, x5, x4, x3");
}

void Task631::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x2, x5, x4, x3") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("x1, c2, x0, x2");
}

void Task632::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, x0, x2") += (*ta1_)("x2, a1") * (*ta2_)("x0, a1, c2, x1") * (-2);
}

void Task633::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x1, x3, x2, x0") * (*ta2_)("x2, x3, x1, x0");
}

void Task634::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, x1, x0") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x2, a1, c2, x3");
}

void Task635::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1, c2, x3") += (*ta1_)("c3, a1, c2, x3") * (*ta2_)("x2, c3") * 2;
}

void Task636::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, x1, x0") += (*ta1_)("c2, x3, x2, a1") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task637::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x0, x1, x4") * (*ta2_)("x1, x0, x5, x4");
}

void Task638::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta1_)("x0, a1, c2, x1") * 2;
}

void Task639::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x0, x1, x2") * (*ta2_)("x1, x0, x3, x2");
}

void Task640::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c3, x2") * (*ta2_)("x1, a1, x0, c3");
}

void Task641::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1, x0, c3") += (*ta1_)("c2, c3") * (*ta2_)("x0, a1, c2, x1") * (-2);
}

void Task642::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a3, c2, x2") * (*ta2_)("x1, c2, x0, a3");
}

void Task643::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, x0, a3") += (*ta1_)("a3, a1") * (*ta2_)("x0, a1, c2, x1") * 2;
}

void Task644::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x2, x3, a1, c2");
}

void Task645::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, c2") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, x2") * 2;
}

void Task646::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta2_)("c2, a1, x0, x1");
}

void Task647::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x0, a1, c2, x1") * e0_ * (-2);
}

void Task648::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x1, a3") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task649::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta2_)("x0, a1, c2, x1");
}

#endif
