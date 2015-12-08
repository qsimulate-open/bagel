//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks8.cc
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

#include <src/smith/caspt2/CASPT2_tasks8.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task350::compute_() {
  (*ta0_)("c2, a3") += (*ta1_)("c2, a3");
}

void Task351::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("c1, x3");
}

void Task352::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2") * 2;
}

void Task353::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("c1, x3");
}

void Task354::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x3") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2") * (-1);
}

void Task355::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a1, x3");
}

void Task356::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, x2") * (-1);
}

void Task357::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a3") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a1, x3");
}

void Task358::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task359::compute_() {
  (*ta0_)("x3, a2") += (*ta1_)("x3, a2");
}

void Task360::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("c1, x5, x3, x4");
}

void Task361::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x3, x4") += (*ta1_)("x5, x3, x2, x4, x1, x0") * (*ta2_)("x0, x1, c1, x2") * (-1);
}

void Task362::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("c1, x5, x4, x3");
}

void Task363::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x0, x1, c1, x2");
}

void Task364::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2") += (*ta1_)("x5, a1, x4, a2") * (*ta2_)("a1, x5, x4, x3");
}

void Task365::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x5, x4, x3") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x0, a1, x1, x2") * 2;
}

void Task366::compute_() {
  (*ta0_)("a2, x1") += (*ta1_)("x1, a2");
}

void Task367::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c1, c3, x1, x0");
}

void Task368::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task369::compute_() {
  (*ta0_)("a2, c1") += (*ta1_)("a2, c1");
}

void Task370::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c3, x0");
}

void Task371::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c3, x1") * 2;
}

void Task372::compute_() {
  (*ta0_)("a2, c3") += (*ta1_)("c3, a2");
}

void Task373::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2") += (*ta1_)("c1, a2, c3, x0") * (*ta2_)("c1, x0");
}

void Task374::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-1);
}

void Task375::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, c1");
}

void Task376::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c1, x0");
}

void Task377::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a4, c1, x0") += (*ta1_)("x1, a4, c1, x0") * 2
     + (*ta1_)("c1, a4, x1, x0") * (-4);
}

void Task378::compute_() {
  (*ta0_)("x1, x2") += (*ta1_)("x2, x1");
}

void Task379::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x2, x1") += (*ta1_)("x0, x3, x2, x1") * (*ta2_)("x0, x3");
}

void Task380::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("c3, a2, c1, x3") * (*ta1_)("c1, a2, c3, x0") * (-1);
}

void Task381::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("c1, a2, c3, x3") * (*ta1_)("c1, a2, c3, x0") * 2;
}

void Task382::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x0, x3");
}

void Task383::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a3, c2, a1") * (*ta1_)("x0, a1, c2, a3") * (-1);
}

void Task384::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, x3") += (*ta1_)("x3, a1, c2, a3") * (*ta1_)("x0, a1, c2, a3") * 2;
}

void Task385::compute_() {
  (*ta0_)("c4, c1") += (*ta1_)("c1, c4");
}

void Task386::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c4") += (*ta1_)("c4, a2, c3, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task387::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task388::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c4") += (*ta1_)("c3, a2, c4, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task389::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0");
}

void Task390::compute_() {
  (*ta0_)("c4, c3") += (*ta1_)("c3, c4");
}

void Task391::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c4") += (*ta1_)("c4, a2, c1, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task392::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0");
}

void Task393::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, c4") += (*ta1_)("c1, a2, c4, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task394::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-2);
}

void Task395::compute_() {
  (*ta0_)("a2, a4") += (*ta1_)("a2, a4");
}

void Task396::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a4") += (*ta1_)("c3, a4, c1, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task397::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * (-1);
}

void Task398::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, a4") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task399::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x0") * 2;
}

#endif
