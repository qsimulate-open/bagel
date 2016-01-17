//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks14.cc
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

#include <src/smith/mrci/MRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task650::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c5, c3, a2") += (*ta1_)("x0, c5, c3, a2")
     + (*ta1_)("x0, a2, c3, c5") * (-2);
}

void Task651::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c5, a2, c1, x1") * (*ta2_)("c5, c3, a4, x1");
}

void Task652::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c5, c3, a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a4");
}

void Task653::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c5, c3, a4") += (*ta1_)("x0, c5, c3, a4") * (-2)
     + (*ta1_)("x0, a4, c3, c5");
}

void Task654::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c5, x1") * (*ta2_)("c5, c3, a2, x1");
}

void Task655::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c5, c3, a2, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a2");
}

void Task656::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c5, c3, a2") += (*ta1_)("x0, c5, c3, a2") * (-2)
     + (*ta1_)("x0, a2, c3, c5");
}

void Task657::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c5, x1") * (*ta2_)("c5, c3, a4, x1");
}

void Task658::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c5, c3, a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a4");
}

void Task659::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x0, c5, c3, a4") += (*ta1_)("x0, c5, c3, a4") * 4
     + (*ta1_)("x0, a4, c3, c5") * (-2);
}

void Task660::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a4, a5, a2") * (*ta2_)("c1, a5, c3, x0");
}

void Task661::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a5, c3, x0") += (*ta1_)("x0, x1") * (*ta2_)("c1, a5, c3, x1") * 2;
}

void Task662::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a2, a5, a4") * (*ta2_)("c1, a5, c3, x0");
}

void Task663::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, a5, c3, x0") += (*ta1_)("x0, x1") * (*ta2_)("c1, a5, c3, x1") * (-1);
}

void Task664::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a4, c1, x2") * (*ta2_)("c3, a2, x3, x2");
}

void Task665::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a2, x1, x0");
}

void Task666::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, x1, x0") += (*ta1_)("c3, a2, x1, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a2") * 0.5;
}

void Task667::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("c3, x1, x0, a2") * 0.5;
}

void Task668::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x1, a2, c3, x0") * (-0.5);
}

void Task669::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a2, c1, x2") * (*ta2_)("c3, a4, x3, x2");
}

void Task670::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a4, x1, x0");
}

void Task671::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a4, x1, x0") += (*ta1_)("c3, a4, x1, x0") * (-1)
     + (*ta1_)("x1, a4, c3, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a4") * (-1);
}

void Task672::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a4") * (-0.5);
}

void Task673::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, c5") * (*ta2_)("a4, c5");
}

void Task674::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, c5") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c5, x0");
}

void Task675::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a4, c5, x0") += (*ta1_)("x1, a4, c5, x0") * 2
     + (*ta1_)("c5, a4, x1, x0") * (-4);
}

void Task676::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c3, c5") * (*ta2_)("a2, c5");
}

void Task677::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c5") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c5, x0");
}

void Task678::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a2, c5, x0") += (*ta1_)("x1, a2, c5, x0") * (-1)
     + (*ta1_)("c5, a2, x1, x0") * 2;
}

void Task679::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a4, a5, a2") * (*ta2_)("a5, c1");
}

void Task680::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a5, c1, x0");
}

void Task681::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a5, c1, x0") += (*ta1_)("x1, a5, c1, x0") * (-2)
     + (*ta1_)("c1, a5, x1, x0") * 4;
}

void Task682::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a2, a5, a4") * (*ta2_)("a5, c1");
}

void Task683::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a5, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a5, c1, x0");
}

void Task684::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a5, c1, x0") += (*ta1_)("x1, a5, c1, x0")
     + (*ta1_)("c1, a5, x1, x0") * (-2);
}

void Task685::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, x3, x2") * (*ta2_)("c3, a2, x3, x2");
}

void Task686::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a2, x1, x0");
}

void Task687::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, x1, x0") += (*ta1_)("c3, a2, x1, x0") * (-1)
     + (*ta1_)("x1, a2, c3, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a2") * (-1);
}

void Task688::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a2") * (-0.5);
}

void Task689::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, x3, x2") * (*ta2_)("c3, a4, x3, x2");
}

void Task690::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a4, x1, x0");
}

void Task691::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a4, x1, x0") += (*ta1_)("c3, a4, x1, x0") * 2
     + (*ta1_)("x1, a4, c3, x0") * (-1)
     + (*ta1_)("x1, x0, c3, a4") * 2;
}

void Task692::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a4");
}

void Task693::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("a4, x0");
}

void Task694::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1");
}

void Task695::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a4") * (*ta2_)("a2, x0");
}

void Task696::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a2, x2, x1") * (-2);
}

void Task697::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a4") * (*ta2_)("c3, c6, a5, a2") * (-8);
}

void Task698::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c6, a5") * (*ta2_)("c3, c6, a5, a2") * 4;
}

void Task699::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a2") * (*ta2_)("c3, c6, a5, a4") * 4;
}

#endif
