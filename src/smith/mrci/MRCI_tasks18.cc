//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks18.cc
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

#include <src/smith/mrci/MRCI_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task850::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a1, x2, x1");
}

void Task851::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a1, x2, x1") += (*ta1_)("a4, a1, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a1") * (-0.5);
}

void Task852::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("a4, x2, x1, a1") * (-0.5);
}

void Task853::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a1, a4, x1") * 0.5;
}

void Task854::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, a4, c2, a1") * (*ta2_)("a4, a3, x3, x0");
}

void Task855::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a3, x2, x1");
}

void Task856::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a3, x2, x1") += (*ta1_)("a4, a3, x2, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a3") * (-0.5);
}

void Task857::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("a4, x2, x1, a3") * (-0.5);
}

void Task858::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a3, a4, x1") * 0.5;
}

void Task859::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, a1, c2, a4") * (*ta2_)("a4, a3, x3, x0");
}

void Task860::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a3, x2, x1");
}

void Task861::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, a3, x2, x1") += (*ta1_)("a4, a3, x2, x1")
     + (*ta1_)("x2, a3, a4, x1") * (-0.5)
     + (*ta1_)("x2, x1, a4, a3");
}

void Task862::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, a3, x3, x0") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("a4, x2, x1, a3") * 0.5;
}

void Task863::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x5, a1, x4, a3") * (*ta2_)("c2, x5, x0, x4");
}

void Task864::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x5, x0, x4") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("c2, x3, x2, x1") * (-1);
}

void Task865::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x5, x0, x4") += (*ta1_)("x5, x0, x4, x1, x3, x2") * (*ta2_)("x3, x2, c2, x1") * (-1);
}

void Task866::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("a4, x1, c2, a1") * (*ta2_)("a3, a4, x0, x1");
}

void Task867::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, a4, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, a4") * (-2);
}

void Task868::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("a4, x1, c2, a3") * (*ta2_)("a1, a4, x0, x1");
}

void Task869::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a4, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a4") * 4;
}

void Task870::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, x1, a4, a1") * (*ta2_)("a3, a4, x1, x0");
}

void Task871::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, a4, x1, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x3, a3, x2, a4") * 2;
}

void Task872::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, x1, a4, a3") * (*ta2_)("a1, a4, x0, x1");
}

void Task873::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a4, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a4") * (-2);
}

void Task874::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, x0") * (*ta2_)("x3, a3, c2, a1");
}

void Task875::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a3, c2, a1") += (*ta1_)("x3, a3, c2, a1")
     + (*ta1_)("x3, a1, c2, a3") * (-2);
}

void Task876::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x5, x0") * (*ta2_)("x5, a3, c2, a1");
}

void Task877::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x5, a3, c2, a1") += (*ta1_)("x5, a3, c2, a1") * 0.5
     + (*ta1_)("x5, a1, c2, a3") * (-1);
}

void Task878::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("a1, x0, x1, a2") + (*ta1_)("a2, x1, x0, a1");
}

void Task879::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x2, a2") * (*ta2_)("a1, x0, x2, x1");
}

void Task880::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
}

void Task881::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x2, x3, a1, a2");
}

void Task882::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a1, c3, a2") * (*ta2_)("x2, c3") * (-1);
}

void Task883::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a1, x2, a3") * (*ta2_)("a3, a2") * 2;
}

void Task884::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a1, c4, a3") * (*ta2_)("x2, c4, a3, a2") * (-1);
}

void Task885::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a4, c3, a1") * (*ta2_)("x2, a2, a4, c3") * (-1);
}

void Task886::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, a2") += (*ta1_)("x3, a1, c3, a4") * (*ta2_)("x2, a2, a4, c3") * 2;
}

void Task887::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x5, a1, c3, x4") * (*ta2_)("a2, c3, x5, x0, x4, x1");
}

void Task888::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x5, x0, x4, x1") += (*ta1_)("x5, x0, x2, x4, x3, x1") * (*ta2_)("x3, a2, x2, c3") * (-1);
}

void Task889::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("c3, a1, x5, x4") * (*ta2_)("a2, c3, x5, x4, x1, x0");
}

void Task890::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x5, x4, x1, x0") += (*ta1_)("x5, x4, x3, x1, x2, x0") * (*ta2_)("x3, a2, x2, c3");
}

void Task891::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x7, a1, x6, x5") * (*ta2_)("a2, x7, x0, x6, x5, x1");
}

void Task892::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x7, x0, x6, x5, x1") += (*ta1_)("x7, x0, x6, x5, x4, x1, x3, x2") * (*ta2_)("x4, a2, x3, x2") * 0.5;
}

void Task893::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x7, x0, x6, x5, x1") += (*ta1_)("x7, x0, x6, x5, x4, x3, x2, x1") * (*ta2_)("x4, x3, x2, a2") * 0.5;
}

void Task894::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x2, a1, a3, a2") * (*ta2_)("a3, x1, x2, x0");
}

void Task895::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x2, x0") += (*ta1_)("x5, x1, x4, x3, x2, x0") * (*ta2_)("x5, a3, x4, x3") * (-1);
}

void Task896::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x5, a1, c3, a2") * (*ta2_)("c3, x5, x0, x1");
}

void Task897::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x0, x1") += (*ta1_)("x5, x0, x4, x1, x3, x2") * (*ta2_)("x4, c3, x3, x2") * (-0.5);
}

void Task898::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x0, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x4, x3, x2, c3") * (-0.5);
}

void Task899::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a2, x3, a1");
}

#endif
