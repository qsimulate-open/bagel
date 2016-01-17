//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks10.cc
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

#include <src/smith/caspt2/CASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task450::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x0, x1") * (*ta2_)("a4, x1");
}

void Task451::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, a4") * 4;
}

void Task452::compute_() {
  (*ta0_)("a4, c3") += (*ta1_)("a4, c3");
}

void Task453::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a2, c1");
}

void Task454::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task455::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-4)
     + (*ta1_)("c1, a2, x1, x0") * 8;
}

void Task456::compute_() {
  (*ta0_)("x0, x1") += (*ta1_)("x1, x0");
}

void Task457::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0") += (*ta1_)("x1, x0") * (*ta2_)("");
}

void Task458::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a4, c3, a2").dot((*ta1_)("c1, a2, c3, a4") * (-4)).get();
}

void Task459::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a2, c3, a4").dot((*ta1_)("c1, a2, c3, a4") * 8).get();
}

void Task460::compute_() {
  (*ta0_)("c5, c3") += (*ta1_)("c3, c5");
}

void Task461::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c1, a4, c5, a2") * (*ta1_)("c1, a2, c3, a4") * 8;
}

void Task462::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c1, a2, c5, a4") * (*ta1_)("c1, a2, c3, a4") * (-16);
}

void Task463::compute_() {
  (*ta0_)("a4, a5") += (*ta1_)("a4, a5");
}

void Task464::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a5") += (*ta1_)("c1, a5, c3, a2") * (*ta1_)("c1, a2, c3, a4") * (-8);
}

void Task465::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a5") += (*ta1_)("c1, a2, c3, a5") * (*ta1_)("c1, a2, c3, a4") * 16;
}

void Task466::compute_() {
  (*ta0_)("x0, c3") += (*ta1_)("c3, x0");
}

void Task467::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x1, x0") * (*ta2_)("c3, x1");
}

void Task468::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, x1") += (*ta1_)("x1, a4, c1, a2") * (*ta1_)("c1, a2, c3, a4") * (-4);
}

void Task469::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, x1") += (*ta1_)("x1, a2, c1, a4") * (*ta1_)("c1, a2, c3, a4") * 2;
}

void Task470::compute_() {
  (*ta0_)("a1, x1") += (*ta1_)("x1, a1");
}

void Task471::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a3, c2, x1, x0");
}

void Task472::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a3, c2, x2");
}

void Task473::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a3, c2, x2") += (*ta1_)("x3, a3, c2, x2") * (-1)
     + (*ta1_)("c2, a3, x3, x2") * 2;
}

void Task474::compute_() {
  (*ta0_)("a3, x1") += (*ta1_)("x1, a3");
}

void Task475::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, c2, x0, x1");
}

void Task476::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2");
}

void Task477::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

void Task478::compute_() {
  (*ta0_)("a1, c2") += (*ta1_)("c2, a1");
}

void Task479::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a3, x0");
}

void Task480::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-1);
}

void Task481::compute_() {
  (*ta0_)("a3, c2") += (*ta1_)("a3, c2");
}

void Task482::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, x0");
}

void Task483::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1") * 2;
}

void Task484::compute_() {
  (*ta0_)("c4, x1") += (*ta1_)("c4, x1");
}

void Task485::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, c4");
}

void Task486::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c4") += (*ta1_)("c2, a3, c4, a1") * (*ta1_)("x0, a1, c2, a3") * (-4);
}

void Task487::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c4") += (*ta1_)("c2, a1, c4, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task488::compute_() {
  (*ta0_)("c4, c2") += (*ta1_)("c2, c4");
}

void Task489::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4") += (*ta1_)("x1, a3, c4, a1") * (*ta2_)("a3, c2, a1, x1");
}

void Task490::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3");
}

void Task491::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4") += (*ta1_)("x1, a1, c4, a3") * (*ta2_)("a3, c2, a1, x1");
}

void Task492::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-2);
}

void Task493::compute_() {
  (*ta0_)("a1, a4") += (*ta1_)("a1, a4");
}

void Task494::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a4") += (*ta1_)("x1, a4, c2, a3") * (*ta2_)("a3, c2, a1, x1");
}

void Task495::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task496::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a4") += (*ta1_)("x1, a3, c2, a4") * (*ta2_)("a3, c2, a1, x1");
}

void Task497::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-1);
}

void Task498::compute_() {
  (*ta0_)("a3, a4") += (*ta1_)("a3, a4");
}

void Task499::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a4") += (*ta1_)("x1, a4, c2, a1") * (*ta2_)("a3, c2, a1, x1");
}

#endif
