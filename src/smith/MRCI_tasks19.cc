//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks19.cc
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

#include <src/smith/MRCI_tasks19.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task900::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x0, x1") += (*ta1_)("x5, x0, x4, x1, x3, x2") * (*ta2_)("x4, c3, x3, x2") * (-0.5);
}

void Task901::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x0, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x4, x3, x2, c3") * (-0.5);
}

void Task902::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, a2, x3, a1");
}

void Task903::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, x3, a1") += (*ta1_)("x3, a3, c4, a1") * (*ta2_)("x2, c4, a3, a2");
}

void Task904::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x5, x0, x4, x1, x3, x2") * (*ta2_)("a2, x3, x2, x5, a1, x4");
}

void Task905::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, x3, x2, x5, a1, x4") += (*ta1_)("x5, a1, x4, a3") * (*ta2_)("a3, a2, x3, x2");
}

void Task906::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a2, x3, x2") += (*ta1_)("a3, a2, x3, x2")
     + (*ta1_)("x3, x2, a3, a2");
}

void Task907::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x3, x2, a2, x5, a1, x4");
}

void Task908::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x2, a2, x5, a1, x4") += (*ta1_)("x5, a1, x4, a3") * (*ta2_)("a3, x3, x2, a2");
}

void Task909::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, a2") += (*ta1_)("x5, x0, x4, x2, x3, x1") * (*ta2_)("x3, a2, x2, x5, a1, x4");
}

void Task910::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, x2, x5, a1, x4") += (*ta1_)("x5, a1, x4, a3") * (*ta2_)("x3, a2, a3, x2") * (-1);
}

void Task911::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x1, x0");
}

void Task912::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, c2, x3, x2");
}

void Task913::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x3, x2") += (*ta1_)("c3, x3, c4, x2") * (*ta2_)("c1, c4, c2, c3") * (-2);
}

void Task914::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x3, x2") += (*ta1_)("c1, a3, c2, a4") * (*ta2_)("a4, x3, a3, x2") * (-2);
}

void Task915::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x0, x5, x1, x4") * (*ta2_)("c1, x5, c2, x4") * (-2);
}

void Task916::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x0, x7, x1, x6") * (*ta2_)("c1, x7, c2, x6") * (-1);
}

void Task917::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, c3, a2, a4");
}

void Task918::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("x1, a2, x0, a4") * (*ta2_)("c1, c3, x1, x0");
}

void Task919::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task920::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c5, a4, c6, a2") * (*ta2_)("c1, c6, c3, c5") * 8;
}

void Task921::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c5, a2, c6, a4") * (*ta2_)("c1, c6, c3, c5") * (-4);
}

void Task922::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c1, a6, c3, a5") * (*ta2_)("a6, a2, a5, a4") * 8;
}

void Task923::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c1, a5, c3, a6") * (*ta2_)("a6, a2, a5, a4") * (-4);
}

void Task924::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("x3, a2, x2, a4") * (*ta2_)("c1, c3, x3, x2");
}

void Task925::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x3, x2") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("c1, x1, c3, x0") * (-2);
}

void Task926::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("") * (*ta2_)("c1, a4, c3, a2");
}

void Task927::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * 4
     + (*ta1_)("c1, a2, c3, a4") * (-8);
}

void Task928::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("") * (*ta2_)("c1, a4, c3, a2");
}

void Task929::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * 2
     + (*ta1_)("c1, a2, c3, a4") * (-4);
}

void Task930::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("x1, x0, a1, a2");
}

void Task931::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("c3, a1, c4, a2") * (*ta2_)("c4, c3, x1, x0");
}

void Task932::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, c3, x1, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x3, c4, x2, c3") * (-2);
}

void Task933::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("a1, a2, x3, x2");
}

void Task934::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a2, x3, x2") += (*ta1_)("x3, a3, x2, a4") * (*ta2_)("a4, a1, a3, a2") * (-2);
}

void Task935::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("x5, x0, x4, x1") * (*ta2_)("x5, a1, x4, a2") * (-2);
}

void Task936::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("x7, x0, x6, x1") * (*ta2_)("x7, a1, x6, a2") * (-1);
}

void Task938::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task939::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
}

void Task940::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3") * 0.5;
}

void Task941::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3") * 0.5;
}

void Task942::compute_() {
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c2, a1, x1, x0");
}

void Task943::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x1, x0") * (*ta2_)("c2, a1") * (-1);
}

void Task944::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2");
}

void Task945::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (-0.5)
     + (*ta1_)("x3, x2, c2, a1") * (-0.5);
}

void Task946::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("c2, x3, x2, a1") * (-0.5);
}

void Task947::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2") * 0.5;
}

void Task948::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("c1, a2, x1, x0");
}

void Task949::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a2, x1, x0") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2") * 2;
}

#endif
