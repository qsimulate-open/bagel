//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks17.cc
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

#include <src/smith/caspt2/CASPT2_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task800::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a1, c2, a3") * (*ta1_)("x0, a1, c2, a3") * 4;
}

void Task801::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x0, x4, x1") * (*ta2_)("x1, x0, x5, x4");
}

void Task802::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a1, x4, a2") * (*ta1_)("x0, a1, x1, a2") * 4;
}

void Task803::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x2, x5, x4, x3, x1, x0") * (*ta2_)("x2, x1, x0, x5, x4, x3");
}

void Task804::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x1, x0, x5, x4, x3") += (*ta1_)("c1, x5, x4, x3") * (*ta2_)("x0, x1, c1, x2");
}

void Task805::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x0, x3, x2, x1") * (*ta2_)("x3, x4, x5, x0, x1, x2");
}

void Task806::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x4, x5, x0, x1, x2") += (*ta1_)("c1, x0, x1, x2") * (*ta2_)("x5, x4, c1, x3");
}

#endif
