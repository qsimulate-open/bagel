//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks10.cc
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

#include <src/smith/caspt2/CASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task450::compute_() {
  (*ta0_)("a2, a3") += (*ta1_)("a2, a3");
}

void Task451::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a3") += (*ta1_)("x3, a3, c1, x2") * (*ta2_)("a2, c1, x3, x2");
}

void Task452::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1") * (-1);
}

void Task453::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a3") += (*ta1_)("c1, a3, x3, x2") * (*ta2_)("a2, c1, x3, x2");
}

void Task454::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1") * 2;
}

void Task455::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a3") += (*ta1_)("x3, a1, x2, a3") * (*ta2_)("a2, a1, x3, x2");
}

void Task456::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, a1, x3, x2") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, a2") * 4;
}

void Task457::compute_() {
  (*ta0_)("x2, c1") += (*ta1_)("x2, c1");
}

void Task458::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c1") += (*ta1_)("c1, a2, x0, x1") * (*ta2_)("a2, x2, x1, x0");
}

void Task459::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a2, x4, x3") * (-1);
}

void Task460::compute_() {
  (*ta0_)("a1, a2") += (*ta1_)("a1, a2");
}

void Task461::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a2") += (*ta1_)("x5, a2, x4, x3") * (*ta2_)("a1, x5, x4, x3");
}

void Task462::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x0, a1, x1, x2");
}

void Task463::compute_() {
  (*ta0_)("a2, x0") += (*ta1_)("a2, x0");
}

void Task464::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x0") += (*ta1_)("x0, x1") * (*ta2_)("a2, x1");
}

void Task465::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta1_)("c1, a2, c3, a4") * (-2);
}

void Task466::compute_() {
  (*ta0_)("a4, x0") += (*ta1_)("a4, x0");
}

void Task467::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x0, x1") * (*ta2_)("a4, x1");
}

void Task468::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x1") += (*ta1_)("c1, a2, c3, x1") * (*ta1_)("c1, a2, c3, a4") * 4;
}

void Task469::compute_() {
  (*ta0_)("a4, c3") += (*ta1_)("a4, c3");
}

void Task470::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a2, c1");
}

void Task471::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task472::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-4)
     + (*ta1_)("c1, a2, x1, x0") * 8;
}

void Task473::compute_() {
  (*ta0_)("x0, x1") += (*ta1_)("x1, x0");
}

void Task474::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0") += (*ta1_)("x1, x0") * (*ta2_)("");
}

void Task475::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a4, c3, a2").dot((*ta1_)("c1, a2, c3, a4") * (-4)).get();
}

void Task476::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("") += (*ta1_)("c1, a2, c3, a4").dot((*ta1_)("c1, a2, c3, a4") * 8).get();
}

void Task477::compute_() {
  (*ta0_)("c5, c3") += (*ta1_)("c3, c5");
}

void Task478::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c1, a4, c5, a2") * (*ta1_)("c1, a2, c3, a4") * 8;
}

void Task479::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c5") += (*ta1_)("c1, a2, c5, a4") * (*ta1_)("c1, a2, c3, a4") * (-16);
}

void Task480::compute_() {
  (*ta0_)("a4, a5") += (*ta1_)("a4, a5");
}

void Task481::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a5") += (*ta1_)("c1, a5, c3, a2") * (*ta1_)("c1, a2, c3, a4") * (-8);
}

void Task482::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a5") += (*ta1_)("c1, a2, c3, a5") * (*ta1_)("c1, a2, c3, a4") * 16;
}

void Task483::compute_() {
  (*ta0_)("x0, c3") += (*ta1_)("c3, x0");
}

void Task484::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x1, x0") * (*ta2_)("c3, x1");
}

void Task485::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, x1") += (*ta1_)("x1, a4, c1, a2") * (*ta1_)("c1, a2, c3, a4") * (-4);
}

void Task486::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, x1") += (*ta1_)("x1, a2, c1, a4") * (*ta1_)("c1, a2, c3, a4") * 2;
}

void Task487::compute_() {
  (*ta0_)("a1, x1") += (*ta1_)("x1, a1");
}

void Task488::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a3, c2, x1, x0");
}

void Task489::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a3, c2, x2");
}

void Task490::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a3, c2, x2") += (*ta1_)("x3, a3, c2, x2") * (-1)
     + (*ta1_)("c2, a3, x3, x2") * 2;
}

void Task491::compute_() {
  (*ta0_)("a3, x1") += (*ta1_)("x1, a3");
}

void Task492::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, c2, x0, x1");
}

void Task493::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2");
}

void Task494::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

void Task495::compute_() {
  (*ta0_)("a1, c2") += (*ta1_)("c2, a1");
}

void Task496::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a3, x0");
}

void Task497::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-1);
}

void Task498::compute_() {
  (*ta0_)("a3, c2") += (*ta1_)("a3, c2");
}

void Task499::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, x0");
}

#endif
