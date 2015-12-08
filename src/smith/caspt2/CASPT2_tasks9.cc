//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks9.cc
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

#include <src/smith/CASPT2_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task400::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c1") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("a2, c3, x1, x0");
}

void Task401::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("x3, a2, c3, x2");
}

void Task402::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, a2, x3, x2") * (-1);
}

void Task403::compute_() {
  (*ta0_)("x1, c3") += (*ta1_)("x1, c3");
}

void Task404::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c3") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("a2, c1, x0, x1");
}

void Task405::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, a2, c1, x2");
}

void Task406::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, c1, x2") += (*ta1_)("x3, a2, c1, x2") * (-1)
     + (*ta1_)("c1, a2, x3, x2") * 2;
}

void Task407::compute_() {
  (*ta0_)("x1, a4") += (*ta1_)("a4, x1");
}

void Task408::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a4");
}

void Task409::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a4, c3, a2") * (*ta1_)("c1, a2, c3, x0") * (-2);
}

void Task410::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a2, c3, a4") * (*ta1_)("c1, a2, c3, x0") * 4;
}

void Task411::compute_() {
  (*ta0_)("a1, x2") += (*ta1_)("x2, a1");
}

void Task412::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("c2, x1, x2, x0");
}

void Task413::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c2, x3") * (-1);
}

void Task414::compute_() {
  (*ta0_)("c3, x2") += (*ta1_)("x2, c3");
}

void Task415::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("c3, a1, c2, x3") * (*ta2_)("c2, a1, x3, x2");
}

void Task416::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x0, a1, c2, x1");
}

void Task417::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("c2, a1, c3, x3") * (*ta2_)("c2, a1, x2, x3");
}

void Task418::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x2, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task419::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("c3, a2, c1, x3") * (*ta2_)("a2, c1, x2, x3");
}

void Task420::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x2, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, a2, x0, x1") * (-1);
}

void Task421::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("c1, a2, c3, x3") * (*ta2_)("a2, c1, x2, x3");
}

void Task422::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x2, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, a2, x0, x1") * 2;
}

void Task423::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("x3, a1, c3, a2") * (*ta2_)("a2, a1, x3, x2");
}

void Task424::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, a1, x3, x2") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, a2") * (-2);
}

void Task425::compute_() {
  (*ta0_)("a1, a3") += (*ta1_)("a1, a3");
}

void Task426::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a3") += (*ta1_)("x3, a3, c2, x2") * (*ta2_)("c2, a1, x3, x2");
}

void Task427::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x0, a1, c2, x1");
}

void Task428::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a3") += (*ta1_)("c2, a3, x3, x2") * (*ta2_)("c2, a1, x3, x2");
}

void Task429::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task430::compute_() {
  (*ta0_)("c3, a4") += (*ta1_)("a4, c3");
}

void Task431::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c2, a4, c3, a1") * (*ta2_)("c2, a1");
}

void Task432::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, x1") * 2;
}

void Task433::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c2, a1, c3, a4") * (*ta2_)("c2, a1");
}

void Task434::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, x1") * (-4);
}

void Task435::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a2, c1");
}

void Task436::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2, x0, x1") * (-4);
}

void Task437::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a2, c1");
}

void Task438::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2, x0, x1") * 8;
}

void Task439::compute_() {
  (*ta0_)("a2, x2") += (*ta1_)("x2, a2");
}

void Task440::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2") += (*ta1_)("c1, a2, x0, x1") * (*ta2_)("c1, x2, x1, x0");
}

void Task441::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3");
}

void Task442::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2") += (*ta1_)("x0, a1, x1, a2") * (*ta2_)("a1, x0, x2, x1");
}

void Task443::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * 2;
}

void Task444::compute_() {
  (*ta0_)("c3, c1") += (*ta1_)("c1, c3");
}

void Task445::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3") += (*ta1_)("x3, a2, c3, x2") * (*ta2_)("a2, c1, x3, x2");
}

void Task446::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1");
}

void Task447::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3") += (*ta1_)("c3, a2, x3, x2") * (*ta2_)("a2, c1, x3, x2");
}

void Task448::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1") * (-2);
}

void Task449::compute_() {
  (*ta0_)("a2, a3") += (*ta1_)("a2, a3");
}

#endif
