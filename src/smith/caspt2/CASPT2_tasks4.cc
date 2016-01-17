//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks4.cc
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

#include <src/smith/caspt2/CASPT2_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task150::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x5, x4, x1, x0") * (*ta2_)("x5, a2, c1, x4");
}

void Task151::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, a2, c1, x4") += (*ta1_)("x5, a2, c1, x4") * (-1)
     + (*ta1_)("c1, a2, x5, x4") * 2;
}

void Task152::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, x3, a2, x2");
}

void Task153::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("x3, a2, c1, x2") * e0_
     + (*ta1_)("c1, a2, x3, x2") * e0_ * (-2);
}

void Task154::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("x3, a2, c3, x2") * (*ta2_)("c1, c3");
}

void Task155::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("x3, a3, c1, x2") * (*ta2_)("a3, a2") * (-1);
}

void Task156::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("c3, a2, x3, x2") * (*ta2_)("c1, c3") * (-2);
}

void Task157::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("c1, a3, x3, x2") * (*ta2_)("a3, a2") * 2;
}

void Task158::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("x3, a3, c1, a2") * (*ta2_)("a3, x2") * 2;
}

void Task159::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3, a2, x2") += (*ta1_)("x3, a2, c1, a3") * (*ta2_)("a3, x2") * (-1);
}

void Task160::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("c1, x2") * (*ta2_)("a2, x2, x1, x0");
}

void Task161::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a2, x4, x3") * (-1);
}

void Task162::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2");
}

void Task163::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a4, c3") * (-4);
}

void Task164::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, c3") * 8;
}

void Task165::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task166::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x3, x4, x2, x1") * (*ta2_)("x3, x5, a1, x4");
}

void Task167::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, a1, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta2_)("x3, c2");
}

void Task168::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a1, x5, x4");
}

void Task169::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1, x5, x4") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("x3, c2") * (-1);
}

void Task170::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x7, x0, x6, x5, x2, x1") * (*ta2_)("x7, a1, x6, x5");
}

void Task171::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("a1, x5, x4, x3");
}

void Task172::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a1, x4, x3") * e0_ * (-1);
}

void Task173::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a2, x4, x3") * (*ta2_)("a2, a1");
}

void Task174::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, a1, x4, a2") * (*ta2_)("a2, x3") * 2;
}

void Task175::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1");
}

void Task176::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a3, c2") * (-1);
}

void Task177::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a1") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, c2") * 2;
}

void Task178::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("a2, c1, a4, c3") + (*ta1_)("a4, c3, a2, c1");
}

void Task179::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("a2, x1");
}

void Task180::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a2") * (-1);
}

void Task181::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, x1") * (*ta2_)("a4, x1");
}

void Task182::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a4") * 2;
}

void Task183::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a2") * (*ta2_)("a4, c1");
}

void Task184::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c1, x0");
}

void Task185::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a4, c1, x0") += (*ta1_)("x1, a4, c1, x0")
     + (*ta1_)("c1, a4, x1, x0") * (-2);
}

void Task186::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a4") * (*ta2_)("a2, c1");
}

void Task187::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task188::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-2)
     + (*ta1_)("c1, a2, x1, x0") * 4;
}

void Task189::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c1, a2") * (*ta2_)("c3, x1");
}

void Task190::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1") += (*ta1_)("x1, x0") * (*ta2_)("c3, x0") * (-2);
}

void Task191::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c1, a4") * (*ta2_)("c3, x1");
}

void Task192::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1") += (*ta1_)("x1, x0") * (*ta2_)("c3, x0");
}

void Task193::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, x0, a1");
}

void Task194::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a1") * (*ta2_)("a3, c2, x1, x0");
}

void Task195::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a3, c2, x2");
}

void Task196::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a3, c2, x2") += (*ta1_)("x3, a3, c2, x2") * (-1)
     + (*ta1_)("c2, a3, x3, x2") * 2;
}

void Task197::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a3") * (*ta2_)("a1, c2, x0, x1");
}

void Task198::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2");
}

void Task199::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

#endif
