//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks3.cc
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

#include <src/smith/relmrci/RelMRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task100::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x5, x3, x4, x1, x2") * (*ta2_)("x3, c2, x2, c1, x5, x4");
}

void Task101::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c2, x2, c1, x5, x4") += (*ta1_)("c1, x5, c3, x4") * (*ta2_)("x3, c3, c2, x2");
}

void Task102::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x7, x6, c1, x5") * (*ta2_)("c2, x7, x6, x0, x5, x1");
}

void Task103::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x6, x0, x5, x1") += (*ta1_)("x7, x6, x0, x5, x1, x4, x3, x2") * (*ta2_)("c2, x4, x3, x2") * 0.5;
}

void Task104::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x6, x0, x5, x1") += (*ta1_)("x7, x6, x0, x5, x4, x3, x1, x2") * (*ta2_)("x4, x3, c2, x2") * 0.5;
}

void Task105::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c1, x2, c2, c3") * (*ta2_)("c3, x1, x0, x2");
}

void Task106::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1, x0, x2") += (*ta1_)("x5, x4, x1, x3, x0, x2") * (*ta2_)("x5, x4, c3, x3");
}

void Task107::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, x1, x5, x0");
}

void Task108::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x5, x0") += (*ta1_)("x1, x5, x0, x4, x3, x2") * (*ta2_)("a3, x4, x3, x2") * (-0.5);
}

void Task109::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x5, x0") += (*ta1_)("x1, x5, x4, x3, x0, x2") * (*ta2_)("x4, x3, a3, x2") * (-0.5);
}

void Task110::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x5, x4, x1, x3, x0, x2") * (*ta2_)("c2, x3, x2, c1, x5, x4");
}

void Task111::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, c1, x5, x4") += (*ta1_)("c1, a3, x5, x4") * (*ta2_)("c2, x3, a3, x2") * (-1);
}

void Task112::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task113::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x3, x4, x1, x0") * (*ta2_)("x3, c1, x5, x4");
}

void Task114::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c1, x5, c2, x4") * (*ta2_)("x3, c2") * 2;
}

void Task115::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c2, x5, c3, x4") * (*ta2_)("x3, c3, c1, c2") * (-2);
}

void Task116::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c2, a3, c1, x5") * (*ta2_)("a3, x4, x3, c2") * 0.5;
}

void Task117::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3");
}

void Task118::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("c1, c2") * (-1);
}

void Task119::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("a2, x3");
}

void Task120::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a3, x3, c1, c2");
}

void Task121::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, x3, c1, c2") += (*ta1_)("a3, x3, c1, c2") * (-1)
     + (*ta1_)("c1, x3, a3, c2");
}

void Task122::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
}

void Task123::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("a3, c2");
}

void Task124::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, c2") * (-1);
}

void Task125::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c4, a3, c2, x3") * (*ta2_)("c1, c4, a3, c2");
}

void Task126::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c4, x3") * (*ta2_)("c1, c4, a3, c2") * (-1);
}

void Task127::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a4, c2, a3") * (*ta2_)("a4, x3, a3, c2") * 2;
}

void Task128::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, a4") * (*ta2_)("a4, x3, a3, c2") * (-2);
}

void Task129::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x7, x5, x6, x4, x3, x1, x0") * (*ta2_)("x5, x4, x3, c1, x7, x6");
}

void Task130::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, c1, x7, x6") += (*ta1_)("c1, x7, c2, x6") * (*ta2_)("x5, c2, x4, x3");
}

void Task131::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x7, x3, x6, x5, x4, x1, x0") * (*ta2_)("x5, x4, x3, c1, x7, x6");
}

void Task132::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, c1, x7, x6") += (*ta1_)("c1, x7, c2, x6") * (*ta2_)("x5, x4, x3, c2");
}

void Task133::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x4, x3, x7, x6, x5");
}

void Task134::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("c1, c2, x4, x3");
}

void Task135::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x4, x3") += (*ta1_)("c1, c2, x4, x3") * (-0.5)
     + (*ta1_)("x4, x3, c1, c2") * (-0.5);
}

void Task136::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("c1, a2, x7, x6") * (*ta2_)("a2, x5, x4, x3") * 0.5;
}

void Task137::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x3, x5, x2, x4, x1, x0") * (*ta2_)("c1, x4, x3, x7, x6, x5");
}

void Task138::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("c1, x4, x3, c2") * (-0.5);
}

void Task139::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x4, x5, x2, x3, x1, x0") * (*ta2_)("x4, c1, x3, x7, x6, x5");
}

void Task140::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, c1, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("x4, c2, c1, x3") * 0.5;
}

void Task141::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

void Task142::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c2, a3, c1, x5") * (*ta2_)("a3, c2, x4, x3");
}

void Task143::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3") * 0.5
     + (*ta1_)("x4, x3, a3, c2") * 0.5;
}

void Task144::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, c2, x4, x3");
}

void Task145::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3") * (-0.5)
     + (*ta1_)("x4, x3, a3, c2") * (-0.5);
}

void Task146::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c3, a2, c1, x5") * (*ta2_)("x4, c3, a2, x3") * (-0.5);
}

void Task147::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x3, x5, x2, x4, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

void Task148::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, x4, x3, c2") * (-0.5);
}

void Task149::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x4, x5, x2, x3, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

#endif
