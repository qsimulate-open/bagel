//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks9.cc
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

#include <src/smith/mrci/MRCI_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task400::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x2, a1, c2, c3") * (*ta2_)("c3, x1, x2, x0");
}

void Task401::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c3, x3");
}

void Task402::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c3, a1, c2, x5") * (*ta2_)("c3, x1, x5, x0");
}

void Task403::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1, x5, x0") += (*ta1_)("x1, x5, x4, x0, x3, x2") * (*ta2_)("x4, c3, x3, x2") * 0.5;
}

void Task404::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1, x5, x0") += (*ta1_)("x1, x5, x4, x3, x2, x0") * (*ta2_)("x4, x3, x2, c3") * 0.5;
}

void Task405::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c2, a1, c3, x5") * (*ta2_)("c3, x5, x1, x0");
}

void Task406::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x1, x0") += (*ta1_)("x4, x5, x3, x2, x1, x0") * (*ta2_)("x4, c3, x3, x2") * (-0.5);
}

void Task407::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("x4, x3, x2, c3") * (-0.5);
}

void Task408::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, a1, c3, x4") * (*ta2_)("c2, c3, x5, x0, x1, x4");
}

void Task409::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x0, x1, x4") += (*ta1_)("x5, x0, x1, x4, x3, x2") * (*ta2_)("c2, c3, x3, x2");
}

void Task410::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3, x3, x2") += (*ta1_)("c2, c3, x3, x2") * (-0.5)
     + (*ta1_)("x3, x2, c2, c3") * (-0.5);
}

void Task411::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x0, x1, x4") += (*ta1_)("x5, x0, x2, x4, x1, x3") * (*ta2_)("c2, x3, x2, c3") * (-0.5);
}

void Task412::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x0, x1, x4") += (*ta1_)("x5, x0, x3, x4, x1, x2") * (*ta2_)("x3, c3, c2, x2") * 0.5;
}

void Task413::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x0, x1, x4, x3, x2") * (*ta2_)("a1, x3, x2, x5, c2, x4");
}

void Task414::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x3, x2, x5, c2, x4") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("a3, a1, x3, x2");
}

void Task415::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a1, x3, x2") += (*ta1_)("a3, a1, x3, x2") * 0.5
     + (*ta1_)("x3, x2, a3, a1") * 0.5;
}

void Task416::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x3, x1, x4, x2, x0") * (*ta2_)("x3, x2, a1, x5, c2, x4");
}

void Task417::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, a1, x5, c2, x4") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("a3, x3, x2, a1") * 0.5;
}

void Task418::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x2, x1, x4, x3, x0") * (*ta2_)("x3, a1, x2, x5, c2, x4");
}

void Task419::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1, x2, x5, c2, x4") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("x3, a1, a3, x2") * (-0.5);
}

void Task420::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c3, a1, x5, x4") * (*ta2_)("c2, c3, x5, x4, x1, x0");
}

void Task421::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x4, x1, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("c2, c3, x3, x2");
}

void Task422::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3, x3, x2") += (*ta1_)("c2, c3, x3, x2") * 0.5
     + (*ta1_)("x3, x2, c2, c3") * 0.5;
}

void Task423::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x4, x1, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("c2, x3, x2, c3") * 0.5;
}

void Task424::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c3, x5, x4, x1, x0") += (*ta1_)("x5, x4, x3, x0, x1, x2") * (*ta2_)("x3, c3, c2, x2") * (-0.5);
}

void Task425::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("a1, x3, x2, c2, x5, x4");
}

void Task426::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x3, x2, c2, x5, x4") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a3, a1, x3, x2");
}

void Task427::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a1, x3, x2") += (*ta1_)("a3, a1, x3, x2") * (-0.5)
     + (*ta1_)("x3, x2, a3, a1") * (-0.5);
}

void Task428::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x3, x2, a1, c2, x5, x4");
}

void Task429::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, a1, c2, x5, x4") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a3, x3, x2, a1") * (-0.5);
}

void Task430::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x4, x3, x0, x1, x2") * (*ta2_)("x3, a1, x2, c2, x5, x4");
}

void Task431::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1, x2, c2, x5, x4") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("x3, a1, a3, x2") * 0.5;
}

void Task432::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x7, a1, x6, x5") * (*ta2_)("c2, x7, x0, x6, x5, x1");
}

void Task433::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x0, x6, x5, x1") += (*ta1_)("x7, x0, x6, x5, x1, x4, x3, x2") * (*ta2_)("c2, x4, x3, x2") * 0.5;
}

void Task434::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x0, x6, x5, x1") += (*ta1_)("x7, x0, x6, x5, x4, x3, x1, x2") * (*ta2_)("x4, x3, c2, x2") * 0.5;
}

void Task435::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("a3, x2, c2, a1") * (*ta2_)("a3, x2, x1, x0");
}

void Task436::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a3, x4, x3") * (-1);
}

void Task437::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c2, x2, a3, a1") * (*ta2_)("a3, x0, x1, x2");
}

void Task438::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0, x1, x2") += (*ta1_)("x5, x0, x4, x3, x1, x2") * (*ta2_)("x5, a3, x4, x3");
}

void Task439::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, a3, c2, a1") * (*ta2_)("a3, x5, x1, x0");
}

void Task440::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x1, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("a3, x4, x3, x2") * (-0.5);
}

void Task441::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x4, x3, a3, x2") * (-0.5);
}

void Task442::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, a1, c2, a3") * (*ta2_)("a3, x5, x0, x1");
}

void Task443::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x0, x1") += (*ta1_)("x5, x0, x1, x4, x3, x2") * (*ta2_)("a3, x4, x3, x2") * 0.5;
}

void Task444::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x5, x0, x1") += (*ta1_)("x5, x0, x4, x3, x1, x2") * (*ta2_)("x4, x3, a3, x2") * 0.5;
}

void Task445::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x0, x4, x2, x1, x3") * (*ta2_)("c2, x3, x2, x5, a1, x4");
}

void Task446::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, x5, a1, x4") += (*ta1_)("x5, a1, x4, a3") * (*ta2_)("c2, x3, a3, x2") * (-2);
}

void Task447::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x0, x1, x4") * (*ta2_)("x5, a1, c2, x4") * (-1);
}

void Task448::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x4, x1, x0") * (*ta2_)("c2, a1, x5, x4");
}

void Task449::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x7, x0, x1, x6") * (*ta2_)("x7, a1, c2, x6") * (-0.5);
}

#endif
