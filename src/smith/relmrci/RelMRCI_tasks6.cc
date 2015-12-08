//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks6.cc
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

#include <src/smith/RelMRCI_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task250::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c4, x0, x3") += (*ta1_)("x2, x3, x0, x1") * (*ta2_)("x2, c4, c3, x1") * (-0.5);
}

void Task251::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, a4, c1, x3") * (*ta2_)("a4, a2, x0, x3");
}

void Task252::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a2, x0, x3") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("a4, a2, x2, x1");
}

void Task253::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a2, x2, x1") += (*ta1_)("a4, a2, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a2") * (-0.5);
}

void Task254::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a2, x0, x3") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("a4, x2, x1, a2") * (-0.5);
}

void Task255::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a2, x0, x3") += (*ta1_)("x2, x3, x0, x1") * (*ta2_)("x2, a2, a4, x1") * 0.5;
}

void Task256::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2, c4, x3") * (*ta2_)("c3, c4, x0, x3");
}

void Task257::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c4, x0, x3") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("c3, c4, x2, x1");
}

void Task258::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c4, x2, x1") += (*ta1_)("c3, c4, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, c3, c4") * (-0.5);
}

void Task259::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c4, x0, x3") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c3, x2, x1, c4") * (-0.5);
}

void Task260::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c4, x0, x3") += (*ta1_)("x2, x3, x0, x1") * (*ta2_)("x2, c4, c3, x1") * 0.5;
}

void Task261::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a4, c3, x3") * (*ta2_)("a4, a2, x0, x3");
}

void Task262::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a2, x0, x3") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("a4, a2, x2, x1");
}

void Task263::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a2, x2, x1") += (*ta1_)("a4, a2, x2, x1") * 0.5
     + (*ta1_)("x2, a2, a4, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a2") * 0.5;
}

void Task264::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a2, x0, x3") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("a4, x2, x1, a2") * 0.5;
}

void Task265::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, a2, x5, x4") * (*ta2_)("c1, x5, x4, x0");
}

void Task266::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x0") += (*ta1_)("x5, x4, x0, x3, x2, x1") * (*ta2_)("c1, x3, x2, x1") * (-0.5);
}

void Task267::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-0.5);
}

void Task268::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("c3, x5, x4, x0");
}

void Task269::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x4, x0") += (*ta1_)("x5, x4, x0, x3, x2, x1") * (*ta2_)("c3, x3, x2, x1") * 0.5;
}

void Task270::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x0, x1") * (*ta2_)("x3, x2, c3, x1") * 0.5;
}

void Task271::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, x1, c1, c4") * (*ta2_)("c4, a2, x0, x1");
}

void Task272::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, a2, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c4, a2, x3, x2") * (-1);
}

void Task273::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("a4, x1, c1, a2") * (*ta2_)("c3, a4, x0, x1");
}

void Task274::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, a4, x3, x2");
}

void Task275::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, x1, c3, c4") * (*ta2_)("c4, a2, x0, x1");
}

void Task276::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, a2, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c4, a2, x3, x2");
}

void Task277::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, x1, a4, a2") * (*ta2_)("c3, a4, x0, x1");
}

void Task278::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, a4, x3, x2") * (-1);
}

void Task279::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("a4, x1, c3, a2") * (*ta2_)("c1, a4, x0, x1");
}

void Task280::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a4, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c1, a4, x3, x2") * (-1);
}

void Task281::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, x1, a4, a2") * (*ta2_)("c1, a4, x0, x1");
}

void Task282::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a4, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c1, a4, x3, x2");
}

void Task283::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, x2, c3, x1") * (*ta2_)("a2, x2, x0, x1");
}

void Task284::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x2, x0, x1") += (*ta1_)("x5, x2, x4, x3, x0, x1") * (*ta2_)("x5, a2, x4, x3") * (-1);
}

void Task285::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x3, a4, c3, a2") * (*ta2_)("c1, a4, x3, x0");
}

void Task286::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a4, x3, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("c1, x2, a4, x1");
}

void Task287::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x3, a2, c3, a4") * (*ta2_)("c1, a4, x3, x0");
}

void Task288::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a4, x3, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c1, x2, a4, x1") * (-1);
}

void Task289::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x3, a4, c1, a2") * (*ta2_)("c3, a4, x3, x0");
}

void Task290::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("c3, x2, a4, x1") * (-1);
}

void Task291::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x3, a2, c1, a4") * (*ta2_)("c3, a4, x3, x0");
}

void Task292::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("c3, x2, a4, x1");
}

void Task293::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x0, x3") * (*ta2_)("c3, a2, c1, x3");
}

void Task294::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x3") += (*ta1_)("c3, a2, c1, x3")
     + (*ta1_)("c1, a2, c3, x3") * (-1);
}

void Task295::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x0, x5") * (*ta2_)("c3, a2, c1, x5");
}

void Task296::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x5") += (*ta1_)("c3, a2, c1, x5") * 0.5
     + (*ta1_)("c1, a2, c3, x5") * (-0.5);
}

void Task297::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("c1, x1, x0, a2");
}

void Task298::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x2, a2") * (*ta2_)("c1, x2, x1, x0");
}

void Task299::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3");
}

#endif
