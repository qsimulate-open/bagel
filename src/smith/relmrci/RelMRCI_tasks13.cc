//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks13.cc
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

#include <src/smith/relmrci/RelMRCI_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task600::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a3, a4, a1") * (*ta2_)("c2, a4, x1, x0");
}

void Task601::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a4, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a4, x3, x2") * (-1);
}

void Task602::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a1, a4, a3") * (*ta2_)("c2, a4, x1, x0");
}

void Task603::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a4, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a4, x3, x2");
}

void Task604::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, x2, x1") * (*ta2_)("a3, x0, x2, x1");
}

void Task605::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
}

void Task606::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, x2, x1") * (*ta2_)("a1, x0, x2, x1");
}

void Task607::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * 0.5;
}

void Task608::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, x2, x1, a1") * (*ta2_)("a3, x2, x1, x0");
}

void Task609::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
}

void Task610::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, x2, x1, a3") * (*ta2_)("a1, x0, x1, x2");
}

void Task611::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, x2") += (*ta1_)("x5, x0, x4, x3, x1, x2") * (*ta2_)("x5, a1, x4, x3") * 0.5;
}

void Task612::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x2, a1, c2, x1") * (*ta2_)("a3, x1, x2, x0");
}

void Task613::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x2, x0") += (*ta1_)("x5, x1, x4, x3, x2, x0") * (*ta2_)("x5, a3, x4, x3") * 0.5;
}

void Task614::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x2, a3, c2, x1") * (*ta2_)("a1, x0, x2, x1");
}

void Task615::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * (-0.5);
}

void Task616::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x2, x1, c2, a1") * (*ta2_)("a3, x0, x2, x1");
}

void Task617::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
}

void Task618::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x2, x1, c2, a3") * (*ta2_)("a1, x0, x2, x1");
}

void Task619::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * 0.5;
}

void Task620::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, a4, a1") * (*ta2_)("a4, x0");
}

void Task621::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1");
}

void Task622::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, a4, a3") * (*ta2_)("a4, x0");
}

void Task623::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1") * (-1);
}

void Task624::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c4, a3, c5, a1") * (*ta2_)("c5, c2, c4, x0");
}

void Task625::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c5, c2, c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c5, c2, c4") * 2;
}

void Task626::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c4, a1, c5, a3") * (*ta2_)("c5, c2, c4, x0");
}

void Task627::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c5, c2, c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c5, c2, c4") * (-2);
}

void Task628::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x3, a3, c4, a1") * (*ta2_)("c2, c4, x3, x0");
}

void Task629::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("c2, c4, x2, x1");
}

void Task630::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4, x2, x1") += (*ta1_)("c2, c4, x2, x1") * 0.5
     + (*ta1_)("x2, x1, c2, c4") * 0.5;
}

void Task631::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, x2, x1, c4") * 0.5;
}

void Task632::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, c4, c2, x1") * (-0.5);
}

void Task633::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x3, a1, c4, a3") * (*ta2_)("c2, c4, x3, x0");
}

void Task634::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("c2, c4, x2, x1");
}

void Task635::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4, x2, x1") += (*ta1_)("c2, c4, x2, x1") * (-0.5)
     + (*ta1_)("x2, c4, c2, x1") * 0.5
     + (*ta1_)("x2, x1, c2, c4") * (-0.5);
}

void Task636::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("c2, x2, x1, c4") * (-0.5);
}

void Task637::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x3, a4, c2, a3") * (*ta2_)("a4, a1, x3, x0");
}

void Task638::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a1, x2, x1");
}

void Task639::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a1, x2, x1") += (*ta1_)("a4, a1, x2, x1") * 0.5
     + (*ta1_)("x2, x1, a4, a1") * 0.5;
}

void Task640::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("a4, x2, x1, a1") * 0.5;
}

void Task641::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a1, a4, x1") * (-0.5);
}

void Task642::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x3, a3, c2, a4") * (*ta2_)("a4, a1, x3, x0");
}

void Task643::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a1, x2, x1");
}

void Task644::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a1, x2, x1") += (*ta1_)("a4, a1, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a1") * (-0.5);
}

void Task645::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("a4, x2, x1, a1") * (-0.5);
}

void Task646::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a1, a4, x1") * 0.5;
}

void Task647::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x3, a4, c2, a1") * (*ta2_)("a4, a3, x3, x0");
}

void Task648::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a3, x2, x1");
}

void Task649::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a3, x2, x1") += (*ta1_)("a4, a3, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a3") * (-0.5);
}

#endif
