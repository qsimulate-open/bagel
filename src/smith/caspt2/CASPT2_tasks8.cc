//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks8.cc
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

#include <src/smith/caspt2/CASPT2_tasks8.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task350::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c1, c3, x1, x0");
}

void Task351::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task352::compute_() {
  (*ta0_)("a2, c1") += (*ta1_)("a2, c1");
}

void Task353::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c3, x0");
}

void Task354::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c3, x1") * 2;
}

void Task355::compute_() {
  (*ta0_)("a2, c3") += (*ta1_)("c3, a2");
}

void Task356::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c1, x0");
}

void Task357::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-1);
}

void Task358::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, c1");
}

void Task359::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c1, x0");
}

void Task360::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a4, c1, x0") += (*ta1_)("x1, a4, c1, x0") * 2
     + (*ta1_)("c1, a4, x1, x0") * (-4);
}

void Task361::compute_() {
  (*ta0_)("x1, x2") += (*ta1_)("x2, x1");
}

void Task362::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x2, x1") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("x0, x3");
}

void Task363::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("c3, a2, c1, x3") * (*ta1_)("c1, a2, c3, x0") * (-1);
}

void Task364::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("c1, a2, c3, x3") * (*ta1_)("c1, a2, c3, x0") * 2;
}

void Task365::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, x3");
}

void Task366::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a3, c2, a1") * (*ta1_)("x0, a1, c2, a3") * (-1);
}

void Task367::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a1, c2, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task368::compute_() {
  (*ta0_)("c4, c1") += (*ta1_)("c1, c4");
}

void Task369::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c4") += (*ta1_)("c4, a2, c3, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task370::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task371::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c4") += (*ta1_)("c3, a2, c4, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task372::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0");
}

void Task373::compute_() {
  (*ta0_)("c4, c3") += (*ta1_)("c3, c4");
}

void Task374::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c4") += (*ta1_)("c4, a2, c1, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task375::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0");
}

void Task376::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c4") += (*ta1_)("c1, a2, c4, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task377::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task378::compute_() {
  (*ta0_)("a2, a4") += (*ta1_)("a2, a4");
}

void Task379::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a4") += (*ta1_)("c3, a4, c1, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task380::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-1);
}

void Task381::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a4") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task382::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * 2;
}

void Task383::compute_() {
  (*ta0_)("x1, c1") += (*ta1_)("x1, c1");
}

void Task384::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c1") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("a2, c3, x1, x0");
}

void Task385::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("x3, a2, c3, x2");
}

void Task386::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, a2, x3, x2") * (-1);
}

void Task387::compute_() {
  (*ta0_)("x1, c3") += (*ta1_)("x1, c3");
}

void Task388::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c3") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("a2, c1, x0, x1");
}

void Task389::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, a2, c1, x2");
}

void Task390::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, c1, x2") += (*ta1_)("x3, a2, c1, x2") * (-1)
     + (*ta1_)("c1, a2, x3, x2") * 2;
}

void Task391::compute_() {
  (*ta0_)("x1, a4") += (*ta1_)("a4, x1");
}

void Task392::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, a4");
}

void Task393::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a4, c3, a2") * (*ta1_)("c1, a2, c3, x0") * (-2);
}

void Task394::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, a4") += (*ta1_)("c1, a2, c3, a4") * (*ta1_)("c1, a2, c3, x0") * 4;
}

void Task395::compute_() {
  (*ta0_)("a1, x2") += (*ta1_)("x2, a1");
}

void Task396::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("c2, x1, x2, x0");
}

void Task397::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c2, x3") * (-1);
}

void Task398::compute_() {
  (*ta0_)("c3, x2") += (*ta1_)("x2, c3");
}

void Task399::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c3") += (*ta1_)("c3, a1, c2, x3") * (*ta2_)("c2, a1, x3, x2");
}

#endif
