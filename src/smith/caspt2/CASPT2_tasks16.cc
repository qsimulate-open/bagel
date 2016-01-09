//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks16.cc
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

#include <src/smith/caspt2/CASPT2_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task750::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0") += (*ta1_)("c2, a1") * (*ta2_)("x1, c2, a1, x0");
}

void Task751::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, a1, x0") += (*ta1_)("x0, a1, c2, x1") * (-2)
     + (*ta1_)("x1, a1, c2, x0") * (-2);
}

void Task752::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0") += (*ta1_)("c1, a2") * (*ta2_)("x1, x0, a2, c1");
}

void Task753::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a2, c1") += (*ta1_)("c1, a2, x0, x1") * 4
     + (*ta1_)("c1, a2, x1, x0") * 4;
}

void Task754::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x2, x4, x3, x1, x0") * (*ta2_)("x1, x0, x2, x5, x4, x3");
}

void Task755::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x2, x5, x4, x3") += (*ta1_)("x5, a2, x4, x3") * (*ta2_)("x1, x0, a2, x2");
}

void Task756::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a2, x2") += (*ta1_)("c1, x2") * (*ta2_)("c1, a2, x0, x1") * (-2);
}

void Task757::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x2, x5, x4, x3") += (*ta1_)("x0, x1, x2, a1") * (*ta2_)("x5, a1, x4, x3");
}

void Task758::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x0, x3, x4, x2, x1") * (*ta2_)("x3, x5, x4, x2, x1, x0");
}

void Task759::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, x4, x2, x1, x0") += (*ta1_)("x0, a1, x1, x2") * (*ta2_)("x3, x5, a1, x4");
}

void Task760::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta2_)("x3, c2") * 2;
}

void Task761::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, x5, x4, x2, x1, x0");
}

void Task762::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, x4, x2, x1, x0") += (*ta1_)("x0, a1, x1, x2") * (*ta2_)("x3, a1, x5, x4");
}

void Task763::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("x3, c2") * (-2);
}

void Task764::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, x4, x2, x1, x0") += (*ta1_)("x5, x4, x3, a1") * (*ta2_)("x0, a1, x1, x2");
}

void Task765::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x7, x0, x6, x5, x2, x1") * (*ta2_)("x2, x1, x0, x7, x6, x5");
}

void Task766::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x7, x6, x5") += (*ta1_)("x7, a1, x6, x5") * (*ta1_)("x0, a1, x1, x2") * 2;
}

void Task767::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x0, x4, x3, x2, x1") * (*ta2_)("x2, x1, x0, x5, x4, x3");
}

void Task768::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a2, x4, x3") * (*ta2_)("x2, x1, x0, a2");
}

void Task769::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, a2") += (*ta1_)("a2, a1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task770::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x0, a1, x1, x2") * (*ta2_)("x3, x5, a1, x4");
}

void Task771::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a1, x4, a2") * (*ta2_)("a2, x3") * 4;
}

void Task772::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * (*ta2_)("x1, a1, x0, x2");
}

void Task773::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1, x0, x2") += (*ta1_)("x0, a1, x1, x2") * e0_ * (-2);
}

void Task774::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1, x0, x2") += (*ta1_)("x2, a2") * (*ta2_)("x0, a1, x1, a2") * 4;
}

void Task775::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * (*ta2_)("x0, a1, x1, x2");
}

void Task776::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x0, a1, x1, x2") * (*ta2_)("x5, a1, x4, x3");
}

void Task777::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x0, x2, x1") * (*ta2_)("x3, x2, x1, x0");
}

void Task778::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x0, a1, x1, x2") * (*ta2_)("x3, a1");
}

void Task779::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a3, c2") * (-2);
}

void Task780::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, c2") * 4;
}

void Task781::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a3, x2, x1") * (*ta2_)("a3, x0");
}

void Task782::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, x0") += (*ta1_)("c2, a1") * (*ta2_)("x0, a1, c2, a3") * (-2);
}

void Task783::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a1, x2, x1") * (*ta2_)("a1, x0");
}

void Task784::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0") += (*ta1_)("c2, a3") * (*ta2_)("x0, a1, c2, a3") * 4;
}

void Task785::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a1, x2, a3") * (*ta2_)("a3, a1, x0, x1");
}

void Task786::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a1, x0, x1") += (*ta1_)("c2, x1") * (*ta2_)("x0, a1, c2, a3") * (-4);
}

void Task787::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a1, x0, x1") += (*ta1_)("a3, a2") * (*ta2_)("x0, a1, x1, a2") * 8;
}

void Task788::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x0, a1, x1, a2") * (*ta2_)("x2, x3, a1, a2");
}

void Task789::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a1, c3, a2") * (*ta2_)("x2, c3") * (-4);
}

void Task790::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a1, x2, a2") * (*ta1_)("x0, a1, x1, a2") * e0_ * (-4);
}

void Task791::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a1, x2, a2") * (*ta2_)("x0, a1, x1, a2") * 2;
}

void Task792::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x0, a1, x1, a2") * (*ta2_)("x3, a1, x2, a2") * 2;
}

void Task793::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x3, a1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task794::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, x1, x0") += (*ta1_)("x0, a1") * (*ta2_)("x3, a1, x2, x1") * 2;
}

void Task795::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0") * (*ta2_)("");
}

void Task796::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a4, c3, a2").dot((*ta1_)("c1, a2, c3, a4") * (-8)).get();
}

void Task797::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a2, c3, a4").dot((*ta1_)("c1, a2, c3, a4") * 16).get();
}

void Task798::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x0") * (*ta2_)("x0, x3");
}

void Task799::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a3, c2, a1") * (*ta1_)("x0, a1, c2, a3") * (-2);
}

#endif
