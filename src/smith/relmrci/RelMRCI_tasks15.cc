//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks15.cc
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

#include <src/smith/relmrci/RelMRCI_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task700::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, c2, x3, x2");
}

void Task701::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x3, x2") += (*ta1_)("c3, x3, c4, x2") * (*ta2_)("c1, c4, c2, c3") * (-2);
}

void Task702::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c2, x3, x2") += (*ta1_)("c1, a3, c2, a4") * (*ta2_)("a4, x3, a3, x2") * (-2);
}

void Task703::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x0, x5, x1, x4") * (*ta2_)("c1, x5, c2, x4") * (-2);
}

void Task704::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x1, x0") += (*ta1_)("x0, x7, x1, x6") * (*ta2_)("c1, x7, c2, x6") * (-1);
}

void Task705::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, c3, a2, a4");
}

void Task706::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("x1, a2, x0, a4") * (*ta2_)("c1, c3, x1, x0");
}

void Task707::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task708::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c5, a4, c6, a2") * (*ta2_)("c1, c6, c3, c5") * 2;
}

void Task709::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("c5, a2, c6, a4") * (*ta2_)("c1, c6, c3, c5") * (-2);
}

void Task710::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("x3, a2, x2, a4") * (*ta2_)("c1, c3, x3, x2");
}

void Task711::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x3, x2") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("c1, x1, c3, x0") * (-2);
}

void Task712::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("") * (*ta2_)("c1, a4, c3, a2");
}

void Task713::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * 2
     + (*ta1_)("c1, a2, c3, a4") * (-2);
}

void Task714::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, a2, a4") += (*ta1_)("") * (*ta2_)("c1, a4, c3, a2");
}

void Task715::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2")
     + (*ta1_)("c1, a2, c3, a4") * (-1);
}

void Task716::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("x1, x0, a1, a2");
}

void Task717::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("c3, a1, c4, a2") * (*ta2_)("c4, c3, x1, x0");
}

void Task718::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c4, c3, x1, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x3, c4, x2, c3") * (-2);
}

void Task719::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("x5, x0, x4, x1") * (*ta2_)("x5, a1, x4, a2") * (-2);
}

void Task720::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("x1, x0, a1, a2") += (*ta1_)("x7, x0, x6, x1") * (*ta2_)("x7, a1, x6, a2") * (-1);
}

void Task722::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task723::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
}

void Task724::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3") * 0.5;
}

void Task725::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("x5, x4, c1, x3") * 0.5;
}

void Task726::compute_() {
  (*ta0_)("x0, x1, c1, a2") += (*ta1_)("c1, a2, x1, x0");
}

void Task727::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a2, x1, x0") += (*ta1_)("x1, x0") * (*ta2_)("c1, a2");
}

void Task728::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c1, a2, x3, x2");
}

void Task729::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, x3, x2") += (*ta1_)("c1, a2, x3, x2") * 0.5
     + (*ta1_)("x3, a2, c1, x2") * (-0.5)
     + (*ta1_)("x3, x2, c1, a2") * 0.5;
}

void Task730::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3, x2, a2") * 0.5;
}

void Task731::compute_() {
  (*ta0_)("x1, x2, x0, a1") += (*ta1_)("a1, x0, x2, x1");
}

void Task732::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1");
}

void Task733::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * 0.5;
}

void Task734::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x5, x4, x3, a1") * 0.5;
}

void Task735::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x0, x1");
}

void Task736::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x3, c2, x2");
}

void Task737::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c3, c1, a2, x0");
}

void Task738::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, c1, a2, x0") += (*ta1_)("x0, x1") * (*ta2_)("c3, x1, c1, a2");
}

void Task739::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, x1, c1, a2") += (*ta1_)("c3, x1, c1, a2")
     + (*ta1_)("c1, x1, c3, a2") * (-1);
}

void Task740::compute_() {
  (*ta0_)("c3, a4, c1, a2") += (*ta1_)("c1, a4, c3, a2");
}

void Task741::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a4, c3, a2") += (*ta1_)("c1, a4, c3, a2") * (-1)
     + (*ta1_)("c1, a2, c3, a4");
}

void Task742::compute_() {
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("a3, c2, a1, x0");
}

void Task743::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, c2, a1, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, a3, c2, a1");
}

void Task744::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a3, c2, a1") += (*ta1_)("x1, a3, c2, a1") * (-1)
     + (*ta1_)("x1, a1, c2, a3");
}

void Task745::compute_() {
  (*ta0_)("x1, a2, x0, a1") += (*ta1_)("a1, a2, x0, x1");
}

void Task746::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a1, a2, x0, x1") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, a2");
}

void Task748::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x0, x1");
}

void Task749::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x3, c2, x2") * 2;
}

#endif
