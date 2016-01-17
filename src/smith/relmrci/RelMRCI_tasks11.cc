//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks11.cc
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

#include <src/smith/relmrci/RelMRCI_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task500::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("a4, x0");
}

void Task501::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1");
}

void Task502::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a4") * (*ta2_)("a2, x0");
}

void Task503::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a2, x2, x1") * (-1);
}

void Task504::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a4") * (*ta2_)("c3, c6, a5, a2") * (-2);
}

void Task505::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c6, a5") * (*ta2_)("c3, c6, a5, a2") * 2;
}

void Task506::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a2") * (*ta2_)("c3, c6, a5, a4") * 2;
}

void Task507::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c6, a5") * (*ta2_)("c3, c6, a5, a4") * (-2);
}

void Task508::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a6, c5, a4") * (*ta2_)("c3, a2, a6, c5") * 2;
}

void Task509::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c5, a6") * (*ta2_)("c3, a2, a6, c5") * (-2);
}

void Task510::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a6, c5, a2") * (*ta2_)("c3, a4, a6, c5") * (-2);
}

void Task511::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c5, a6") * (*ta2_)("c3, a4, a6, c5") * 2;
}

void Task512::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a4, c1, a2") * (*ta2_)("c3, x3");
}

void Task513::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, x2, x1, x0") * (-0.5);
}

void Task514::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x2, x1, c3, x0") * (-0.5);
}

void Task515::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a2, c1, a4") * (*ta2_)("c3, x3");
}

void Task516::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, x2, x1, x0") * 0.5;
}

void Task517::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x2, x1, c3, x0") * 0.5;
}

void Task518::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c5, a2") * (*ta2_)("c1, c3, c5, x1");
}

void Task519::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, c5, x1") += (*ta1_)("x1, x0") * (*ta2_)("c1, x0, c3, c5") * (-1);
}

void Task520::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c5, a4") * (*ta2_)("c1, c3, c5, x1");
}

void Task521::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, c5, x1") += (*ta1_)("x1, x0") * (*ta2_)("c1, x0, c3, c5");
}

void Task522::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a5, c1, a4") * (*ta2_)("a5, c3, a2, x1");
}

void Task523::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c3, a2, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a2");
}

void Task524::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, x0, c3, a2") += (*ta1_)("a5, x0, c3, a2") * (-1)
     + (*ta1_)("c3, x0, a5, a2");
}

void Task525::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c1, a5") * (*ta2_)("a5, c3, a2, x1");
}

void Task526::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c3, a2, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a2");
}

void Task527::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, x0, c3, a2") += (*ta1_)("a5, x0, c3, a2")
     + (*ta1_)("c3, x0, a5, a2") * (-1);
}

void Task528::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a5, c1, a2") * (*ta2_)("a5, c3, a4, x1");
}

void Task529::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c3, a4, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a4");
}

void Task530::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, x0, c3, a4") += (*ta1_)("a5, x0, c3, a4")
     + (*ta1_)("c3, x0, a5, a4") * (-1);
}

void Task531::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c1, a5") * (*ta2_)("a5, c3, a4, x1");
}

void Task532::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c3, a4, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a4");
}

void Task533::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a5, x0, c3, a4") += (*ta1_)("a5, x0, c3, a4") * (-1)
     + (*ta1_)("c3, x0, a5, a4");
}

void Task534::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, x0, a1");
}

void Task535::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a1") * (*ta2_)("c2, a3, x1, x0");
}

void Task536::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a3, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a3, x3, x2");
}

void Task537::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a3") * (*ta2_)("c2, a1, x1, x0");
}

void Task538::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

void Task539::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1") * (*ta2_)("a3, x0");
}

void Task540::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-1);
}

void Task541::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3") * (*ta2_)("a1, x0");
}

void Task542::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1");
}

void Task543::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, c4, a1") * (*ta2_)("c4, x0");
}

void Task544::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c4") * (-2);
}

void Task545::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, c4, x2, x1") * (-1);
}

void Task546::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, c4") * (-1);
}

void Task547::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, c4, a3") * (*ta2_)("c4, x0");
}

void Task548::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c4") * 2;
}

void Task549::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, c4, x2, x1");
}

#endif
