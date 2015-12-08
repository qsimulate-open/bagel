//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks14.cc
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

#include <src/smith/CASPT2_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task650::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, x1, x0") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x2, a1, c2, x3");
}

void Task651::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1, c2, x3") += (*ta1_)("c3, a1, c2, x3") * (*ta2_)("x2, c3") * 2;
}

void Task652::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, x1, x0") += (*ta1_)("c2, x3, x2, a1") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task653::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x0, x1, x4") * (*ta2_)("x1, x0, x5, x4");
}

void Task654::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a1, c2, x4") * (*ta1_)("x0, a1, c2, x1") * 2;
}

void Task655::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x0, x1, x2") * (*ta2_)("x1, x0, x3, x2");
}

void Task656::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c3, x2") * (*ta2_)("x1, a1, x0, c3");
}

void Task657::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1, x0, c3") += (*ta1_)("c2, c3") * (*ta2_)("x0, a1, c2, x1") * (-2);
}

void Task658::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a3, c2, x2") * (*ta2_)("x1, c2, x0, a3");
}

void Task659::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, x0, a3") += (*ta1_)("a3, a1") * (*ta2_)("x0, a1, c2, x1") * 2;
}

void Task660::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x2, x3, a1, c2");
}

void Task661::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, a1, c2") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, x2") * 2;
}

void Task662::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta2_)("c2, a1, x0, x1");
}

void Task663::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x0, a1, c2, x1") * e0_ * (-2);
}

void Task664::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x1, a3") * (*ta2_)("x0, a1, c2, a3") * 2;
}

void Task665::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a1, c2, x2") * (*ta2_)("x0, a1, c2, x1");
}

void Task666::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x3, a1, c2, x2");
}

void Task667::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x5, x4, x1, x0") * (*ta2_)("x1, x0, x5, x4");
}

void Task668::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c2, a1, x5, x4") * (*ta1_)("x0, a1, c2, x1") * (-2);
}

void Task669::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta1_)("c1, a2, x0, x1") * (-2);
}

void Task670::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x5, x4") += (*ta1_)("c1, a2, x5, x4") * (*ta1_)("c1, a2, x0, x1") * 4;
}

void Task671::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("ci0") += (*ta1_)("ci0, x3, x2, x1, x0") * (*ta2_)("x1, x0, x3, x2");
}

void Task672::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c3, a1, x3, x2") * (*ta2_)("x1, a1, x0, c3");
}

void Task673::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, a1, x0, c3") += (*ta1_)("c2, c3") * (*ta2_)("x0, a1, c2, x1") * 2;
}

void Task674::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a3, x3, x2") * (*ta2_)("x1, c2, x0, a3");
}

void Task675::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, x0, a3") += (*ta1_)("a3, a1") * (*ta2_)("x0, a1, c2, x1") * (-2);
}

void Task676::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, c2, x0, a3") += (*ta1_)("x1, a1") * (*ta2_)("x0, a1, c2, a3") * 4;
}

void Task677::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x0, a1, c2, x1") * (*ta2_)("x2, x3, c2, a1");
}

void Task678::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, c2, a1") += (*ta1_)("x3, a3, c2, a1") * (*ta2_)("a3, x2") * (-2);
}

void Task679::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a2, c3, x2") * (*ta2_)("x1, x0, a2, c3");
}

void Task680::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a2, c3") += (*ta1_)("c1, c3") * (*ta2_)("c1, a2, x0, x1") * 2;
}

void Task681::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a3, c1, x2") * (*ta2_)("x1, x0, c1, a3");
}

void Task682::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, c1, a3") += (*ta1_)("a3, a2") * (*ta2_)("c1, a2, x0, x1") * (-2);
}

void Task683::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c3, a2, x3, x2") * (*ta2_)("x1, x0, a2, c3");
}

void Task684::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, a2, c3") += (*ta1_)("c1, c3") * (*ta2_)("c1, a2, x0, x1") * (-4);
}

void Task685::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a3, x3, x2") * (*ta2_)("x1, x0, c1, a3");
}

void Task686::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, c1, a3") += (*ta1_)("a3, a2") * (*ta2_)("c1, a2, x0, x1") * 4;
}

void Task687::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x0, x1") * (*ta2_)("x2, x3, c1, a2");
}

void Task688::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, c1, a2") += (*ta1_)("x3, a3, c1, a2") * (*ta2_)("a3, x2") * 4;
}

void Task689::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, x3, c1, a2") += (*ta1_)("x3, a2, c1, a3") * (*ta2_)("a3, x2") * (-2);
}

void Task690::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a3, c2, x2") * (*ta2_)("a3, c2, x0, x1");
}

void Task691::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, c2, x0, x1") += (*ta1_)("x1, a1") * (*ta2_)("x0, a1, c2, a3") * (-2);
}

void Task692::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (*ta2_)("c2, a1, x0, x1");
}

void Task693::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x0, a1, c2, x1") * e0_ * 2;
}

void Task694::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, a1, x0, x1") += (*ta1_)("x1, a3") * (*ta2_)("x0, a1, c2, a3") * (-2);
}

void Task695::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, a2, c1, x2") * (*ta1_)("c1, a2, x0, x1") * e0_ * 2;
}

void Task696::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x3, x2") * (*ta1_)("c1, a2, x0, x1") * e0_ * (-4);
}

void Task697::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c2, a1, x3, x2") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task698::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("x3, x2, c2, a1") * (*ta2_)("x0, a1, c2, x1") * (-1);
}

void Task699::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x1, x0, x3, x2") += (*ta1_)("c1, a2, x3, x2") * (*ta2_)("c1, a2, x0, x1") * 2;
}

#endif
