//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks11.cc
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

#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task500::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * (-1);
}

void Task501::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, a4") += (*ta1_)("x1, a1, c2, a4") * (*ta2_)("a3, c2, a1, x1");
}

void Task502::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x1") += (*ta1_)("x1, x0") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task503::compute_() {
  (*ta0_)("x1, c2") += (*ta1_)("x1, c2");
}

void Task504::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2") += (*ta1_)("x0, a1, c2, a3") * (*ta2_)("a1, a3, x0, x1");
}

void Task505::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a3, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a3") * (-2);
}

void Task507::compute_() {
  (*ta0_)("c1, x0") += (*ta1_)("c1, x0");
}

void Task508::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1");
}

void Task509::compute_() {
  (*ta0_)("c1, a2") += (*ta1_)("a2, c1");
}

void Task510::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c1, x0");
}

void Task511::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c1, x0") += (*ta1_)("x1, a2, c1, x0") * (-1)
     + (*ta1_)("c1, a2, x1, x0") * 2;
}

void Task512::compute_() {
  (*ta0_)("x0, a1") += (*ta1_)("a1, x0");
}

void Task513::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1");
}

void Task515::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x0, x1");
}

void Task516::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x3, c2, x2") * 2;
}

void Task517::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task518::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3");
}

void Task519::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c3, a2, c1, x0");
}

void Task520::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, c1, x0") += (*ta1_)("x0, x1") * (*ta2_)("c3, a2, c1, x1");
}

void Task521::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x1") += (*ta1_)("c3, a2, c1, x1") * (-1)
     + (*ta1_)("c1, a2, c3, x1") * 2;
}

void Task522::compute_() {
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("a1, c2, x0, x1");
}

void Task523::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c2, x2");
}

void Task524::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
}

void Task525::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("a2, c1, x1, x0");
}

void Task526::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a2, c1, x2");
}

void Task527::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, c1, x2") += (*ta1_)("x3, a2, c1, x2") * (-1)
     + (*ta1_)("c1, a2, x3, x2") * 2;
}

void Task528::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task529::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
}

void Task530::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, a4, c3, a2");
}

void Task531::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * (-4)
     + (*ta1_)("c1, a2, c3, a4") * 8;
}

void Task532::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, a1, x0");
}

void Task533::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, a3, c2, a1");
}

void Task534::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3, c2, a1") += (*ta1_)("x1, a3, c2, a1") * (-1)
     + (*ta1_)("x1, a1, c2, a3") * 2;
}

void Task535::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("a1, a2, x0, x1");
}

void Task536::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a2, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a2") * 2;
}

void Task538::compute_() {
  (*ta0_)("ci0") += (*ta1_)("ci0");
}

void Task539::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x0, x5, x1, x4") * (*ta2_)("x1, x0, x5, x4");
}

void Task540::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c1, x5, c2, x4") * (*ta1_)("c1, x0, c2, x1") * 4;
}

void Task541::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x0, x3, x1, x2") * (*ta2_)("x1, x0, x3, x2");
}

void Task542::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, c3, x2") * (*ta2_)("x1, x0, c1, c3");
}

void Task543::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, c1, c3") += (*ta1_)("c2, c3") * (*ta2_)("c1, x0, c2, x1") * (-8);
}

void Task544::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, c2, x2") * (*ta1_)("c1, x0, c2, x1") * e0_ * (-4);
}

void Task545::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x3, c2, x2") * (*ta2_)("c1, x0, c2, x1") * 2;
}

void Task546::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, x0, c2, x1") * (*ta2_)("c1, x3, c2, x2") * 2;
}

void Task547::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x0, x3, x1, x2") * (*ta2_)("x1, x0, x2, x5, x4, x3");
}

void Task548::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x2, x5, x4, x3") += (*ta1_)("x5, x4, c1, x3") * (*ta2_)("x1, x0, c1, x2");
}

void Task549::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, c1, x2") += (*ta1_)("c2, x2") * (*ta2_)("c1, x0, c2, x1") * 4;
}

#endif
