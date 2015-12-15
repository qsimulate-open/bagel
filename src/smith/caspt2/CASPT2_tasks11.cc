//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks11.cc
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

#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task500::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1") * 2;
}

void Task501::compute_() {
  (*ta0_)("c4, x1") += (*ta1_)("c4, x1");
}

void Task502::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, c4");
}

void Task503::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c4") += (*ta1_)("c2, a3, c4, a1") * (*ta1_)("x0, a1, c2, a3") * (-4);
}

void Task504::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c4") += (*ta1_)("c2, a1, c4, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task505::compute_() {
  (*ta0_)("c4, c2") += (*ta1_)("c2, c4");
}

void Task506::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4") += (*ta1_)("x1, a3, c4, a1") * (*ta2_)("a3, c2, a1, x1");
}

void Task507::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3");
}

void Task508::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c4") += (*ta1_)("x1, a1, c4, a3") * (*ta2_)("a3, c2, a1, x1");
}

void Task509::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-2);
}

void Task510::compute_() {
  (*ta0_)("a1, a4") += (*ta1_)("a1, a4");
}

void Task511::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a4") += (*ta1_)("x1, a4, c2, a3") * (*ta2_)("a3, c2, a1, x1");
}

void Task512::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task513::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a1, a4") += (*ta1_)("x1, a3, c2, a4") * (*ta2_)("a3, c2, a1, x1");
}

void Task514::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-1);
}

void Task515::compute_() {
  (*ta0_)("a3, a4") += (*ta1_)("a3, a4");
}

void Task516::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a4") += (*ta1_)("x1, a4, c2, a1") * (*ta2_)("a3, c2, a1, x1");
}

void Task517::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-1);
}

void Task518::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a4") += (*ta1_)("x1, a1, c2, a4") * (*ta2_)("a3, c2, a1, x1");
}

void Task519::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task520::compute_() {
  (*ta0_)("x1, c2") += (*ta1_)("x1, c2");
}

void Task521::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, a3, x0, x1");
}

void Task522::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a3, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a3") * (-2);
}

void Task524::compute_() {
  (*ta0_)("c1, x0") += (*ta1_)("c1, x0");
}

void Task525::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1");
}

void Task526::compute_() {
  (*ta0_)("c1, a2") += (*ta1_)("a2, c1");
}

void Task527::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task528::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-1)
     + (*ta1_)("c1, a2, x1, x0") * 2;
}

void Task529::compute_() {
  (*ta0_)("x0, a1") += (*ta1_)("a1, x0");
}

void Task530::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1");
}

void Task532::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x0, x1");
}

void Task533::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x3, c2, x2") * 2;
}

void Task534::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task535::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3");
}

void Task536::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c3, a2, c1, x0");
}

void Task537::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x0") += (*ta1_)("x0, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task538::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("c3, a2, c1, x1") * (-1)
     + (*ta1_)("c1, a2, c3, x1") * 2;
}

void Task539::compute_() {
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("a1, c2, x0, x1");
}

void Task540::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2");
}

void Task541::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

void Task542::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("a2, c1, x1, x0");
}

void Task543::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a2, c1, x2");
}

void Task544::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, c1, x2") += (*ta1_)("x3, a2, c1, x2") * (-1)
     + (*ta1_)("c1, a2, x3, x2") * 2;
}

void Task545::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task546::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
}

void Task547::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, a4, c3, a2");
}

void Task548::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * (-4)
     + (*ta1_)("c1, a2, c3, a4") * 8;
}

void Task549::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, a1, x0");
}

#endif
