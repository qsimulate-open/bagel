//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks7.cc
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

#include <src/smith/relmrci/RelMRCI_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task300::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c3, a2, c1, x3") * (*ta2_)("x2, c3") * (-1);
}

void Task301::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c1, a2, c3, x3") * (*ta2_)("x2, c3");
}

void Task302::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c4, a2, c3, x3") * (*ta2_)("x2, c4, c1, c3");
}

void Task303::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c3, a2, c4, x3") * (*ta2_)("x2, c4, c1, c3") * (-1);
}

void Task304::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c4, a3, c1, x3") * (*ta2_)("x2, c4, a3, a2") * (-1);
}

void Task305::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c3, a4, c1, x3") * (*ta2_)("x2, a2, a4, c3");
}

void Task306::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c1, a3, c4, x3") * (*ta2_)("x2, c4, a3, a2");
}

void Task307::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c1, a4, c3, x3") * (*ta2_)("x2, a2, a4, c3") * (-1);
}

void Task308::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a4, x3, x2, c3") * (-1);
}

void Task309::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a2, c1, x3") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, x3, x2, c3");
}

void Task310::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x3, x2");
}

void Task311::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c3, a2, x3, x2") * (*ta2_)("c1, c3") * (-1);
}

void Task312::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a3, x3, x2") * (*ta2_)("a3, a2");
}

void Task313::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a3, c1, a2") * (*ta2_)("a3, x2");
}

void Task314::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a2, c1, a3") * (*ta2_)("a3, x2") * (-1);
}

void Task315::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c4, a3, x3, x2") * (*ta2_)("c1, c4, a3, a2") * (-1);
}

void Task316::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c3, a4, x3, x2") * (*ta2_)("c1, a2, a4, c3");
}

void Task317::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a4, c3, x3, x2");
}

void Task318::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3, x3, x2") += (*ta1_)("a4, c3, x3, x2") * (-1)
     + (*ta1_)("x3, x2, a4, c3") * (-1);
}

void Task319::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, c3, x3, x2");
}

void Task320::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, c3, x3, x2") += (*ta1_)("a4, c3, x3, x2")
     + (*ta1_)("x3, x2, a4, c3");
}

void Task321::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a3, c4, a2") * (*ta2_)("x3, c4, a3, x2");
}

void Task322::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a2, c4, a3") * (*ta2_)("x3, c4, a3, x2") * (-1);
}

void Task323::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a4, c3, a2") * (*ta2_)("a4, x2, c1, c3");
}

void Task324::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x2, c1, c3") += (*ta1_)("a4, x2, c1, c3") * (-1)
     + (*ta1_)("c1, x2, a4, c3");
}

void Task325::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a2, c3, a4") * (*ta2_)("a4, x2, c1, c3");
}

void Task326::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x2, c1, c3") += (*ta1_)("a4, x2, c1, c3")
     + (*ta1_)("c1, x2, a4, c3") * (-1);
}

void Task327::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a4, c1, a3") * (*ta2_)("a4, x2, a3, a2");
}

void Task328::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("x3, a3, c1, a4") * (*ta2_)("a4, x2, a3, a2") * (-1);
}

void Task329::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("c1, x2") * (*ta2_)("a2, x2, x1, x0");
}

void Task330::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a2, x4, x3") * (-1);
}

void Task331::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2");
}

void Task332::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a4, c3") * (-2);
}

void Task333::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, c3") * 2;
}

void Task334::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c3, a4, c5, a2") * (*ta2_)("c1, c5, a4, c3") * (-2);
}

void Task335::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c3, a2, c5, a4") * (*ta2_)("c1, c5, a4, c3") * 2;
}

void Task336::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a5, c3, a4") * (*ta2_)("a5, a2, a4, c3") * 2;
}

void Task337::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2") += (*ta1_)("c1, a4, c3, a5") * (*ta2_)("a5, a2, a4, c3") * (-2);
}

void Task338::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x3, a2, x2, c3") * (*ta2_)("c1, c3, x3, x2, x1, x0");
}

void Task339::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x3, x2, x1, x0") += (*ta1_)("x3, x5, x2, x4, x1, x0") * (*ta2_)("c1, x5, c3, x4") * (-2);
}

void Task340::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x4, a2, x3, x2") * (*ta2_)("c1, x4, x3, x2, x1, x0");
}

void Task341::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x4, x3, x2, x1, x0") += (*ta1_)("x7, x6, x4, x5, x3, x2, x1, x0") * (*ta2_)("x7, x6, c1, x5") * 0.5;
}

void Task342::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x4, x3, x2, a2") * (*ta2_)("c1, x2, x4, x3, x1, x0");
}

void Task343::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x4, x3, x1, x0") += (*ta1_)("x7, x6, x2, x5, x4, x3, x1, x0") * (*ta2_)("x7, x6, c1, x5") * 0.5;
}

void Task344::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x2, c3, c1, a2") * (*ta2_)("c3, x2, x1, x0");
}

void Task345::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c3, x3");
}

void Task346::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("x2, a2, c1, c3") * (*ta2_)("c3, x2, x1, x0");
}

void Task347::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c3, x3") * (-1);
}

void Task348::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x1, x0, a2") += (*ta1_)("c3, a2, c1, x5") * (*ta2_)("c3, x5, x1, x0");
}

void Task349::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x5, x1, x0") += (*ta1_)("x4, x5, x3, x2, x1, x0") * (*ta2_)("x4, c3, x3, x2") * (-0.5);
}

#endif
