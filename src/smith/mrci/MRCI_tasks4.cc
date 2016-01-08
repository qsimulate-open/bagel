//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks4.cc
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

#include <src/smith/mrci/MRCI_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task150::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
}

void Task151::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("a3, c2") * 2;
}

void Task152::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, c2") * (-1);
}

void Task153::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c4, a3, c2, x3") * (*ta2_)("c1, c4, a3, c2");
}

void Task154::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c4, x3") * (*ta2_)("c1, c4, a3, c2") * (-2);
}

void Task155::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a4, c2, a3") * (*ta2_)("a4, x3, a3, c2") * 4;
}

void Task156::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, a4") * (*ta2_)("a4, x3, a3, c2") * (-2);
}

void Task157::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x3, x2, x4, x1, x0") * (*ta2_)("x3, x5, c1, x4");
}

void Task158::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("a2, x3") * (-1);
}

void Task159::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("a3, x3, c1, c2");
}

void Task160::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a2, c1, a3") * (*ta2_)("a3, x4, a2, x3") * (-1);
}

void Task161::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x7, x5, x6, x4, x3, x1, x0") * (*ta2_)("x5, x4, x3, c1, x7, x6");
}

void Task162::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, c1, x7, x6") += (*ta1_)("c1, x7, c2, x6") * (*ta2_)("x5, c2, x4, x3");
}

void Task163::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x7, x3, x6, x5, x4, x1, x0") * (*ta2_)("x5, x4, x3, c1, x7, x6");
}

void Task164::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, c1, x7, x6") += (*ta1_)("c1, x7, c2, x6") * (*ta2_)("x5, x4, x3, c2");
}

void Task165::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x4, x3, x7, x6, x5");
}

void Task166::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("c1, c2, x4, x3");
}

void Task167::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x4, x3") += (*ta1_)("c1, c2, x4, x3") * (-0.5)
     + (*ta1_)("x4, x3, c1, c2") * (-0.5);
}

void Task168::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("c1, a2, x7, x6") * (*ta2_)("a2, x5, x4, x3") * 0.5;
}

void Task169::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x3, x5, x2, x4, x1, x0") * (*ta2_)("c1, x4, x3, x7, x6, x5");
}

void Task170::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("c1, x4, x3, c2") * (-0.5);
}

void Task171::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x4, x5, x2, x3, x1, x0") * (*ta2_)("x4, c1, x3, x7, x6, x5");
}

void Task172::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, c1, x3, x7, x6, x5") += (*ta1_)("x7, x6, c2, x5") * (*ta2_)("x4, c2, c1, x3") * 0.5;
}

void Task173::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

void Task174::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c2, a3, c1, x5") * (*ta2_)("a3, c2, x4, x3");
}

void Task175::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3")
     + (*ta1_)("x4, x3, a3, c2");
}

void Task176::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, c2, x4, x3");
}

void Task177::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x4, x3") += (*ta1_)("a3, c2, x4, x3") * (-0.5)
     + (*ta1_)("x4, x3, a3, c2") * (-0.5);
}

void Task178::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c3, a2, c1, x5") * (*ta2_)("x4, c3, a2, x3") * (-0.5);
}

void Task179::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x3, x5, x2, x4, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

void Task180::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, x4, x3, c2") * (-0.5);
}

void Task181::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x4, x5, x2, x3, x1, x0") * (*ta2_)("x4, x3, c1, x5");
}

void Task182::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x4, x3, c1, x5") += (*ta1_)("c1, a2, c3, x5") * (*ta2_)("x4, c3, a2, x3") * 0.5;
}

void Task183::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x5, x2, x6, x4, x3, x1, x0") * (*ta2_)("x5, x4, x3, x7, c1, x6");
}

void Task184::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, x7, c1, x6") += (*ta1_)("x7, a2, c1, x6") * (*ta2_)("a2, x5, x4, x3") * (-0.5);
}

void Task185::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x3, x2, x6, x5, x4, x1, x0") * (*ta2_)("x5, x4, x3, x7, c1, x6");
}

void Task186::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, x7, c1, x6") += (*ta1_)("x7, a2, c1, x6") * (*ta2_)("x5, x4, a2, x3") * (-0.5);
}

void Task187::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, x3, c1, x7, x6");
}

void Task188::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, x4, x3, c1, x7, x6") += (*ta1_)("c1, a2, x7, x6") * (*ta2_)("x5, x4, a2, x3") * 0.5;
}

void Task189::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x3, x6, x5, x2, x4, x1, x0") * (*ta2_)("c1, x4, x3, x7, x6, x5");
}

void Task190::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x4, x3, x7, x6, x5") += (*ta1_)("x7, a2, x6, x5") * (*ta2_)("c1, x4, a2, x3") * (-1);
}

void Task191::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x2, x5, x1, x0") * (*ta2_)("x7, x6, c1, x5") * (-1);
}

void Task192::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x9, x8, x2, x7, x1, x0") * (*ta2_)("x9, x8, c1, x7") * (-0.5);
}

void Task193::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c1, c3, x0, a2");
}

void Task194::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x1, a2") * (*ta2_)("c1, c3, x1, x0");
}

void Task195::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task196::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2") * (*ta2_)("c3, x0");
}

void Task197::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c3, x1") * 2;
}

void Task198::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, a2") * (*ta2_)("c1, x0");
}

void Task199::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-1);
}

#endif
