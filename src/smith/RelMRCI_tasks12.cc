//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks12.cc
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

#include <src/smith/RelMRCI_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task550::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c4") * 2;
}

void Task551::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, c4, x2, x1");
}

void Task552::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, c4");
}

void Task553::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, x0") * (*ta2_)("c2, x1, a3, a1");
}

void Task554::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c4, a1") * (*ta2_)("c2, c4");
}

void Task555::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c4, a3") * (*ta2_)("c2, c4") * (-1);
}

void Task556::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a3") * (*ta2_)("a4, a1");
}

void Task557::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c2, a4") * (*ta2_)("a4, a1") * (-1);
}

void Task558::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a1") * (*ta2_)("a4, a3") * (-1);
}

void Task559::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c2, a4") * (*ta2_)("a4, a3");
}

void Task560::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a4, c5, a3") * (*ta2_)("x1, c5, a4, a1") * 2;
}

void Task561::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a3, c5, a4") * (*ta2_)("x1, c5, a4, a1") * (-2);
}

void Task562::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a4, c5, a1") * (*ta2_)("x1, c5, a4, a3") * (-2);
}

void Task563::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a1, c5, a4") * (*ta2_)("x1, c5, a4, a3") * 2;
}

void Task564::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a5, c4, a3") * (*ta2_)("x1, a1, a5, c4") * (-2);
}

void Task565::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a3, c4, a5") * (*ta2_)("x1, a1, a5, c4") * 2;
}

void Task566::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a5, c4, a1") * (*ta2_)("x1, a3, a5, c4") * 2;
}

void Task567::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a1, c4, a5") * (*ta2_)("x1, a3, a5, c4") * (-2);
}

void Task568::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c5, a3") * (*ta2_)("c2, c5, a4, a1") * (-1);
}

void Task569::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c5, a4") * (*ta2_)("c2, c5, a4, a1");
}

void Task570::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c5, a1") * (*ta2_)("c2, c5, a4, a3");
}

void Task571::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c5, a4") * (*ta2_)("c2, c5, a4, a3") * (-1);
}

void Task572::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c4, a3") * (*ta2_)("c2, a1, a5, c4");
}

void Task573::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c4, a5") * (*ta2_)("c2, a1, a5, c4") * (-1);
}

void Task574::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c4, a1") * (*ta2_)("c2, a3, a5, c4") * (-1);
}

void Task575::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c4, a5") * (*ta2_)("c2, a3, a5, c4");
}

void Task576::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c2, a4") * (*ta2_)("a5, a1, a4, a3");
}

void Task577::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a5") * (*ta2_)("a5, a1, a4, a3") * (-1);
}

void Task578::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, x1") * (*ta2_)("a1, a3, x0, x1");
}

void Task579::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a3, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a3") * (-2);
}

void Task580::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x2, a1, x1, a3") * (*ta2_)("c2, x1, x2, x0");
}

void Task581::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c2, x3") * (-1);
}

void Task582::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c4, a3, c2, x3") * (*ta2_)("a1, c4, x3, x0");
}

void Task583::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a1, x1, c4");
}

void Task584::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c4, a1, c2, x3") * (*ta2_)("a3, c4, x3, x0");
}

void Task585::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x3, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x2, a3, x1, c4") * (-1);
}

void Task586::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, c4, x3") * (*ta2_)("a1, c4, x3, x0");
}

void Task587::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a1, x1, c4") * (-1);
}

void Task588::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, c4, x3") * (*ta2_)("a3, c4, x3, x0");
}

void Task589::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a3, x1, c4");
}

void Task590::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a1, x5, x4, x0");
}

void Task591::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1") * 0.5;
}

void Task592::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, a1") * 0.5;
}

void Task593::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("a3, x5, x4, x0");
}

void Task594::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x4, x0") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-0.5);
}

void Task595::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, a3") * (-0.5);
}

void Task596::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, c4, c2, a1") * (*ta2_)("c4, a3, x1, x0");
}

void Task597::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, a3, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c4, a3, x3, x2");
}

void Task598::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, c4, c2, a3") * (*ta2_)("c4, a1, x1, x0");
}

void Task599::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, a1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c4, a1, x3, x2") * (-1);
}

#endif
