//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks16.cc
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

#include <src/smith/MRCI_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task750::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a3") * (*ta2_)("a4, a1") * 2;
}

void Task751::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c2, a4") * (*ta2_)("a4, a1") * (-1);
}

void Task752::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a1") * (*ta2_)("a4, a3") * (-1);
}

void Task753::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c2, a4") * (*ta2_)("a4, a3") * 2;
}

void Task754::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a4, c5, a3") * (*ta2_)("x1, c5, a4, a1") * 2;
}

void Task755::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a3, c5, a4") * (*ta2_)("x1, c5, a4, a1") * (-4);
}

void Task756::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a4, c5, a1") * (*ta2_)("x1, c5, a4, a3") * (-4);
}

void Task757::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a1, c5, a4") * (*ta2_)("x1, c5, a4, a3") * 2;
}

void Task758::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a5, c4, a3") * (*ta2_)("x1, a1, a5, c4") * (-4);
}

void Task759::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a3, c4, a5") * (*ta2_)("x1, a1, a5, c4") * 8;
}

void Task760::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a5, c4, a1") * (*ta2_)("x1, a3, a5, c4") * 2;
}

void Task761::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("c2, a1, c4, a5") * (*ta2_)("x1, a3, a5, c4") * (-4);
}

void Task762::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c5, a3") * (*ta2_)("c2, c5, a4, a1") * (-2);
}

void Task763::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c5, a4") * (*ta2_)("c2, c5, a4, a1");
}

void Task764::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c5, a1") * (*ta2_)("c2, c5, a4, a3");
}

void Task765::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c5, a4") * (*ta2_)("c2, c5, a4, a3") * (-2);
}

void Task766::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c4, a3") * (*ta2_)("c2, a1, a5, c4");
}

void Task767::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a3, c4, a5") * (*ta2_)("c2, a1, a5, c4") * (-2);
}

void Task768::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c4, a1") * (*ta2_)("c2, a3, a5, c4") * (-2);
}

void Task769::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a1, c4, a5") * (*ta2_)("c2, a3, a5, c4") * 4;
}

void Task770::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a5, c2, a4") * (*ta2_)("a5, a1, a4, a3") * 2;
}

void Task771::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, a3, a1") += (*ta1_)("x1, a4, c2, a5") * (*ta2_)("a5, a1, a4, a3") * (-1);
}

void Task772::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, x1") * (*ta2_)("a1, a3, x0, x1");
}

void Task773::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a3, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a3") * (-2);
}

void Task774::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x2, a1, x1, a3") * (*ta2_)("c2, x1, x2, x0");
}

void Task775::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c2, x3") * (-1);
}

void Task776::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c4, a3, c2, x3") * (*ta2_)("a1, c4, x3, x0");
}

void Task777::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a1, x1, c4");
}

void Task778::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c4, a1, c2, x3") * (*ta2_)("a3, c4, x3, x0");
}

void Task779::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x3, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x2, a3, x1, c4") * (-1);
}

void Task780::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a3, c4, x3") * (*ta2_)("a1, c4, x3, x0");
}

void Task781::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a1, x1, c4") * (-2);
}

void Task782::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a1, c4, x3") * (*ta2_)("a3, c4, x3, x0");
}

void Task783::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x3, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a3, x1, c4");
}

void Task784::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("a1, x5, x4, x0");
}

void Task785::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1") * (-0.5);
}

void Task786::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, a1") * (-0.5);
}

void Task787::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x5, a1, c2, x4") * (*ta2_)("a3, x5, x0, x4");
}

void Task788::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x0, x4") += (*ta1_)("x5, x0, x3, x4, x2, x1") * (*ta2_)("x3, a3, x2, x1") * 0.5;
}

void Task789::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x0, x4") += (*ta1_)("x5, x0, x1, x4, x3, x2") * (*ta2_)("x3, x2, x1, a3") * 0.5;
}

void Task790::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, c4, c2, a1") * (*ta2_)("a3, c4, x0, x1");
}

void Task791::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a3, c4, x2") * (-1);
}

void Task792::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c4, a3, x3, x2");
}

void Task793::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, c4, c2, a3") * (*ta2_)("a1, c4, x0, x1");
}

void Task794::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c4, x2") * 2;
}

void Task795::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c4, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c4, a1, x3, x2") * (-2);
}

void Task796::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a1, c2, c4") * (*ta2_)("a3, c4, x1, x0");
}

void Task797::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c4, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a3, c4, x2");
}

void Task798::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a3, c4, x2") += (*ta1_)("x3, a3, c4, x2")
     + (*ta1_)("c4, a3, x3, x2") * (-2);
}

void Task799::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a3, c2, c4") * (*ta2_)("a1, c4, x0, x1");
}

#endif
