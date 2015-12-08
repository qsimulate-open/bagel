//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks6.cc
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

#include <src/smith/CASPT2_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task250::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (*ta2_)("x0, a1, c2, x1") * (-0.5);
}

void Task251::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, x2, c2, a1") * (*ta2_)("x0, a1, c2, x1") * (-0.5);
}

void Task252::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x3, x2") * (*ta2_)("c1, a2, x0, x1");
}

void Task253::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a2, c1, x2") * (*ta2_)("c1, a2, x0, x1") * (-0.5);
}

void Task254::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, x2, c1, a2") * (*ta2_)("c1, a2, x0, x1");
}

void Task255::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x1, x3, x2, x0").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task256::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, x3, x2, a1") * (*ta2_)("x0, a1, c2, x1") * (-0.5);
}

void Task257::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x1, x2").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task258::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta2_)("x0, a1, c2, x1") * 0.5;
}

void Task259::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x2, x3, x1, x0").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task260::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, x2, a2") * (*ta2_)("c1, a2, x0, x1") * 0.5;
}

void Task261::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x5, x0, x4, x3, x2, x1").dot((*ta2_)("x2, x1, x0, x5, x4, x3")).get();
}

void Task262::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * (*ta2_)("x0, a1, x1, x2") * 0.5;
}

void Task263::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x5, x4, x3, x0, x2, x1").dot((*ta2_)("x2, x1, x0, x5, x4, x3")).get();
}

void Task264::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, x4, x3, a1") * (*ta2_)("x0, a1, x1, x2") * 0.5;
}

void Task265::compute_() {
  target_ += (*ta1_)("c1, a4, c3, a2").dot((*ta2_)("c1, a2, c3, a4") * (-2)).get();
}

void Task266::compute_() {
  target_ += (*ta1_)("c1, a2, c3, a4").dot((*ta2_)("c1, a2, c3, a4") * 4).get();
}

void Task267::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x1, x0").dot((*ta2_)("x0, x1")).get();
}

void Task268::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a3, c2, a1") * (*ta2_)("x0, a1, c2, a3") * (-1);
}

void Task269::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a1, c2, a3") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task270::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c2, a1") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task271::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2") * (*ta2_)("c1, a2, x0, x1") * 2;
}

void Task272::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x2, x1").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task273::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, x2, a2") * (*ta2_)("x0, a1, x1, a2");
}

void Task274::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1") * (*ta2_)("x0, a1, x1, x2");
}

void Task275::compute_() {
  target_ += (*ta1_)("c1, x3").dot((*ta2_)("c1, x3")).get();
}

void Task276::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2");
}

void Task277::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x0, x3, x1, x2").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task278::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, c2, x2") * (*ta1_)("c1, x0, c2, x1") * 2;
}

void Task279::compute_() {
  target_ += (*ta1_)("x5, x4, c1, x3").dot((*ta2_)("c1, x5, x4, x3")).get();
}

void Task280::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2");
}

void Task281::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x0, x1").dot((*ta2_)("x0, x1")).get();
}

void Task282::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c1, x1") * (*ta1_)("c1, a2, c3, x0") * (-1);
}

void Task283::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, x0") * 2;
}

void Task284::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x1, x2").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task285::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta1_)("x0, a1, c2, x1");
}

void Task286::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x2, x1, x0").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task287::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (*ta1_)("x0, a1, c2, x1") * (-1);
}

void Task288::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a2, c1, x2") * (*ta1_)("c1, a2, x0, x1") * (-1);
}

void Task289::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x3, x2") * (*ta1_)("c1, a2, x0, x1") * 2;
}

void Task290::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x5, x0, x4, x3, x2, x1").dot((*ta2_)("x2, x1, x0, x5, x4, x3")).get();
}

void Task291::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * (*ta1_)("x0, a1, x1, x2");
}

void Task292::compute_() {
  target_ += (*ta1_)("c1, a4, c3, a2").dot((*ta1_)("c1, a2, c3, a4") * (-4)).get();
}

void Task293::compute_() {
  target_ += (*ta1_)("c1, a2, c3, a4").dot((*ta1_)("c1, a2, c3, a4") * 8).get();
}

void Task294::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x1, x0").dot((*ta2_)("x0, x1")).get();
}

void Task295::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a3, c2, a1") * (*ta1_)("x0, a1, c2, a3") * (-1);
}

void Task296::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a1, c2, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task297::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x2, x1").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task298::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, x2, a2") * (*ta1_)("x0, a1, x1, a2") * 2;
}

#endif
