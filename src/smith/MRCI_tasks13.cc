//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks13.cc
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

#include <src/smith/MRCI_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task600::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x0, x6, x3, x5, x4, x2, x1") * (*ta2_)("x5, x4, x3, x7, a1, x6");
}

void Task601::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, x7, a1, x6") += (*ta1_)("x7, a1, x6, a2") * (*ta2_)("x5, x4, a2, x3");
}

void Task602::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x0, x6, x5, x2, x1") * (*ta2_)("x7, a1, x6, x5") * (-1);
}

void Task603::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x9, x0, x8, x7, x2, x1") * (*ta2_)("x9, a1, x8, x7") * (-0.5);
}

void Task604::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("a2, c1, a4, c3") + (*ta1_)("a4, c3, a2, c1");
}

void Task605::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("a2, x1");
}

void Task606::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a2") * (-1);
}

void Task607::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, x1") * (*ta2_)("a4, x1");
}

void Task608::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a4") * 2;
}

void Task609::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a2") * (*ta2_)("a4, c1");
}

void Task610::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c1, x0");
}

void Task611::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a4, c1, x0") += (*ta1_)("x1, a4, c1, x0")
     + (*ta1_)("c1, a4, x1, x0") * (-2);
}

void Task612::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a4") * (*ta2_)("a2, c1");
}

void Task613::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task614::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-2)
     + (*ta1_)("c1, a2, x1, x0") * 4;
}

void Task615::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c5, a2") * (*ta2_)("c3, c5");
}

void Task616::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c3, c5") * 4;
}

void Task617::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c5") += (*ta1_)("x1, x0") * (*ta2_)("c3, c5, x1, x0");
}

void Task618::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5, x1, x0") += (*ta1_)("c3, c5, x1, x0") * 2
     + (*ta1_)("x1, c5, c3, x0") * (-1)
     + (*ta1_)("x1, x0, c3, c5") * 2;
}

void Task619::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c5") += (*ta1_)("x0, x1") * (*ta2_)("c3, x1, x0, c5");
}

void Task620::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c5, a4") * (*ta2_)("c3, c5");
}

void Task621::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c3, c5") * (-8);
}

void Task622::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c5") += (*ta1_)("x1, x0") * (*ta2_)("c3, c5, x1, x0");
}

void Task623::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5, x1, x0") += (*ta1_)("c3, c5, x1, x0") * (-4)
     + (*ta1_)("x1, c5, c3, x0") * 2
     + (*ta1_)("x1, x0, c3, c5") * (-4);
}

void Task624::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c5") += (*ta1_)("x0, x1") * (*ta2_)("c3, x1, x0, c5") * (-2);
}

void Task625::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c3, a2") * (*ta2_)("a5, a4");
}

void Task626::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, a4") += (*ta1_)("a5, a4") * (-4);
}

void Task627::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, a4") += (*ta1_)("x1, x0") * (*ta2_)("a5, a4, x1, x0");
}

void Task628::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, a4, x1, x0") += (*ta1_)("a5, a4, x1, x0") * (-2)
     + (*ta1_)("x1, a4, a5, x0")
     + (*ta1_)("x1, x0, a5, a4") * (-2);
}

void Task629::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, a4") += (*ta1_)("x0, x1") * (*ta2_)("a5, x1, x0, a4") * (-1);
}

void Task630::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, a5") * (*ta2_)("a5, a4");
}

void Task631::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, a4") += (*ta1_)("a5, a4") * 8;
}

void Task632::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, a4") += (*ta1_)("x1, x0") * (*ta2_)("a5, a4, x1, x0");
}

void Task633::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, a4, x1, x0") += (*ta1_)("a5, a4, x1, x0") * 4
     + (*ta1_)("x1, a4, a5, x0") * (-2)
     + (*ta1_)("x1, x0, a5, a4") * 4;
}

void Task634::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, a4") += (*ta1_)("x0, x1") * (*ta2_)("a5, x1, x0, a4") * 2;
}

void Task635::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c1, a2") * (*ta2_)("c3, x1");
}

void Task636::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1") += (*ta1_)("x1, x0") * (*ta2_)("c3, x0") * (-2);
}

void Task637::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c1, a4") * (*ta2_)("c3, x1");
}

void Task638::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1") += (*ta1_)("x1, x0") * (*ta2_)("c3, x0");
}

void Task639::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a4, c3, a2") * (*ta2_)("c1, x0");
}

void Task640::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-1);
}

void Task641::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a2, c3, a4") * (*ta2_)("c1, x0");
}

void Task642::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * 2;
}

void Task643::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c3, x3") * (*ta2_)("a2, x3");
}

void Task644::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x2, a2, x1, x0") * (-0.5);
}

void Task645::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x3") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("x2, x1, x0, a2") * (-0.5);
}

void Task646::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, x3") * (*ta2_)("a4, x3");
}

void Task647::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x2, a4, x1, x0");
}

void Task648::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x3") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("x2, x1, x0, a4");
}

void Task649::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c5, a4, c1, x1") * (*ta2_)("c5, c3, a2, x1");
}

#endif
