//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks16.cc
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

#include <src/smith/relmrci/RelMRCI_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task750::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task751::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3");
}

void Task752::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c3, a2, c1, x0");
}

void Task753::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x0") += (*ta1_)("x0, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task754::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("c3, a2, c1, x1") * (-1)
     + (*ta1_)("c1, a2, c3, x1");
}

void Task755::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("c1, a2, x1, x0");
}

void Task756::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x3, x2");
}

void Task757::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task758::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
}

void Task759::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, a4, c3, a2");
}

void Task760::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * (-2)
     + (*ta1_)("c1, a2, c3, a4") * 2;
}

void Task761::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, a1, x0");
}

void Task762::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, a3, c2, a1");
}

void Task763::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3, c2, a1") += (*ta1_)("x1, a3, c2, a1") * (-1)
     + (*ta1_)("x1, a1, c2, a3");
}

void Task764::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("a1, a2, x0, x1");
}

void Task765::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a2, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a2") * 2;
}

#endif
