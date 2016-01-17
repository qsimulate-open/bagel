//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks6.cc
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

#include <src/smith/caspt2/CASPT2_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task250::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task251::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * 0.5;
}

void Task252::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x5, x4, x3, a1") * 0.5;
}

void Task253::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1");
}

void Task254::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, a4, c3, a2");
}

void Task255::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * (-2)
     + (*ta1_)("c1, a2, c3, a4") * 4;
}

void Task256::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, a1, x0");
}

void Task257::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, a3, c2, a1");
}

void Task258::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3, c2, a1") += (*ta1_)("x1, a3, c2, a1") * (-1)
     + (*ta1_)("x1, a1, c2, a3") * 2;
}

void Task259::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("a1, a2, x0, x1");
}

void Task260::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a2, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a2");
}

void Task261::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x0, x3, x1, x2").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task262::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, c2, x2") * (*ta1_)("c1, x0, c2, x1") * 2;
}

void Task263::compute_() {
  target_ += (*ta1_)("x5, x4, c1, x3").dot((*ta2_)("c1, x5, x4, x3")).get();
}

void Task264::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2");
}

void Task265::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x0, x1").dot((*ta2_)("x0, x1")).get();
}

void Task266::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c3, a2, c1, x1") * (*ta1_)("c1, a2, c3, x0") * (-1);
}

void Task267::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, x0") * 2;
}

void Task268::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x1, x2").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task269::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta1_)("x0, a1, c2, x1");
}

void Task270::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x2, x1, x0").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task271::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (*ta1_)("x0, a1, c2, x1") * (-1);
}

void Task272::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a2, c1, x2") * (*ta1_)("c1, a2, x0, x1") * (-1);
}

void Task273::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x3, x2") * (*ta1_)("c1, a2, x0, x1") * 2;
}

void Task274::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x5, x0, x4, x3, x2, x1").dot((*ta2_)("x2, x1, x0, x5, x4, x3")).get();
}

void Task275::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * (*ta1_)("x0, a1, x1, x2");
}

void Task276::compute_() {
  target_ += (*ta1_)("c1, a4, c3, a2").dot((*ta1_)("c1, a2, c3, a4") * (-4)).get();
}

void Task277::compute_() {
  target_ += (*ta1_)("c1, a2, c3, a4").dot((*ta1_)("c1, a2, c3, a4") * 8).get();
}

void Task278::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x1, x0").dot((*ta2_)("x0, x1")).get();
}

void Task279::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a3, c2, a1") * (*ta1_)("x0, a1, c2, a3") * (-1);
}

void Task280::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x1") += (*ta1_)("x1, a1, c2, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task281::compute_() {
  ta1_->init();
  target_ += (*ta1_)("x3, x0, x2, x1").dot((*ta2_)("x1, x0, x3, x2")).get();
}

void Task282::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, x2, a2") * (*ta1_)("x0, a1, x1, a2") * 2;
}

void Task284::compute_() {
  (*ta0_)("x2, x3") += (*ta1_)("x3, x2");
}

void Task285::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x3, x2") += (*ta1_)("x0, x5, x1, x4, x3, x2") * (*ta2_)("x1, x0, x5, x4");
}

void Task286::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c1, x5, c2, x4") * (*ta1_)("c1, x0, c2, x1") * 2;
}

void Task287::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x3, x2") += (*ta1_)("x5, x0, x1, x4, x3, x2") * (*ta2_)("x1, x0, x5, x4");
}

void Task288::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta1_)("x0, a1, c2, x1");
}

void Task289::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x3, x2") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x1, x0, x5, x4");
}

void Task290::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c2, a1, x5, x4") * (*ta1_)("x0, a1, c2, x1") * (-1);
}

void Task291::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta1_)("c1, a2, x0, x1") * (-1);
}

void Task292::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c1, a2, x5, x4") * (*ta1_)("c1, a2, x0, x1") * 2;
}

void Task293::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x3, x2") += (*ta1_)("x5, x0, x4, x1, x3, x2") * (*ta2_)("x1, x0, x5, x4");
}

void Task294::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a1, x4, a2") * (*ta1_)("x0, a1, x1, a2") * 2;
}

void Task295::compute_() {
  (*ta0_)("c3, c2") += (*ta1_)("c2, c3");
}

void Task296::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3") += (*ta1_)("c1, x3, c3, x2") * (*ta2_)("c2, c1, x3, x2");
}

void Task297::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x0, c2, x1") * (-4);
}

void Task298::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3") += (*ta1_)("x3, a1, c3, x2") * (*ta2_)("c2, a1, x3, x2");
}

void Task299::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

#endif
