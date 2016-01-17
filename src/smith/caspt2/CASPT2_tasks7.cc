//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks7.cc
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

#include <src/smith/caspt2/CASPT2_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task300::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3") += (*ta1_)("c3, a1, x3, x2") * (*ta2_)("c2, a1, x3, x2");
}

void Task301::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x0, a1, c2, x1");
}

void Task302::compute_() {
  (*ta0_)("x2, c2") += (*ta1_)("x2, c2");
}

void Task303::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c2") += (*ta1_)("c1, x0, c2, x1") * (*ta2_)("c1, x0, x1, x2");
}

void Task304::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0, x1, x2") += (*ta1_)("x5, x4, x0, x3, x1, x2") * (*ta2_)("x5, x4, c1, x3") * 2;
}

void Task305::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c2") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("a1, x0, x1, x2");
}

void Task306::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x1, x2") += (*ta1_)("x5, x0, x4, x3, x1, x2") * (*ta2_)("x5, a1, x4, x3");
}

void Task307::compute_() {
  (*ta0_)("x2, a3") += (*ta1_)("x2, a3");
}

void Task308::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("c2, c1, x3, x2");
}

void Task309::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x0, c2, x1") * (-2);
}

void Task310::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a3") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("c2, a1, x3, x2");
}

void Task311::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task312::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a3") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("c2, a1, x3, x2");
}

void Task313::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, a1, x3, x2") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x0, a1, c2, x1");
}

void Task314::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a3") += (*ta1_)("x3, a3, c1, a2") * (*ta2_)("a2, c1, x3, x2");
}

void Task315::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1") * 2;
}

void Task316::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a3") += (*ta1_)("x3, a2, c1, a3") * (*ta2_)("a2, c1, x3, x2");
}

void Task317::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x0, x1") * (-1);
}

void Task318::compute_() {
  (*ta0_)("c2, x3") += (*ta1_)("x3, c2");
}

void Task319::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c2") += (*ta1_)("c1, x5, c2, x4") * (*ta2_)("c1, x5, x3, x4");
}

void Task320::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x3, x4") += (*ta1_)("x2, x5, x3, x4, x1, x0") * (*ta2_)("x0, x1, c1, x2") * 2;
}

void Task321::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c2") += (*ta1_)("x5, a1, c2, x4") * (*ta2_)("a1, x5, x3, x4");
}

void Task322::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x3, x4") += (*ta1_)("x5, x0, x3, x4, x2, x1") * (*ta2_)("x0, a1, x1, x2");
}

void Task323::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c2") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("a1, x5, x4, x3");
}

void Task324::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, x2") * (-1);
}

void Task325::compute_() {
  (*ta0_)("x3, x4") += (*ta1_)("x4, x3");
}

void Task326::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x4, x3") += (*ta1_)("x7, x6, x2, x5, x4, x3, x1, x0") * (*ta2_)("x2, x1, x0, x7, x6, x5");
}

void Task327::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x7, x6, x5") += (*ta1_)("x7, x6, c1, x5") * (*ta1_)("x0, x1, c1, x2");
}

void Task328::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x4, x3") += (*ta1_)("x7, x0, x6, x5, x4, x3, x2, x1") * (*ta2_)("x2, x1, x0, x7, x6, x5");
}

void Task329::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x7, x6, x5") += (*ta1_)("x7, a1, x6, x5") * (*ta1_)("x0, a1, x1, x2");
}

void Task330::compute_() {
  (*ta0_)("c2, c1") += (*ta1_)("c1, c2");
}

void Task331::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("c1, x5, x4, x3");
}

void Task332::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2") * (-1);
}

void Task333::compute_() {
  (*ta0_)("c2, a3") += (*ta1_)("c2, a3");
}

void Task334::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("c1, x3");
}

void Task335::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2") * 2;
}

void Task336::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("c1, x3");
}

void Task337::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2") * (-1);
}

void Task338::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a1, x3");
}

void Task339::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, x2") * (-1);
}

void Task340::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a1, x3");
}

void Task341::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task342::compute_() {
  (*ta0_)("x3, a2") += (*ta1_)("x3, a2");
}

void Task343::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("c1, x5, x3, x4");
}

void Task344::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x3, x4") += (*ta1_)("x5, x3, x2, x4, x1, x0") * (*ta2_)("x0, x1, c1, x2") * (-1);
}

void Task345::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("c1, x5, x4, x3");
}

void Task346::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2");
}

void Task347::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("x5, a1, x4, a2") * (*ta2_)("a1, x5, x4, x3");
}

void Task348::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task349::compute_() {
  (*ta0_)("a2, x1") += (*ta1_)("x1, a2");
}

#endif
