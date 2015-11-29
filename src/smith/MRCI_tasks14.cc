//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks14.cc
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

#include <src/smith/MRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task650::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c3, a2, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a2");
  madness::World::get_default().gop.fence();
}

void Task651::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c5, c3, a2") += (*ta1_)("x0, c5, c3, a2")
     + (*ta1_)("x0, a2, c3, c5") * (-2);
  madness::World::get_default().gop.fence();
}

void Task652::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c5, a2, c1, x1") * (*ta2_)("c5, c3, a4, x1");
  madness::World::get_default().gop.fence();
}

void Task653::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c3, a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a4");
  madness::World::get_default().gop.fence();
}

void Task654::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c5, c3, a4") += (*ta1_)("x0, c5, c3, a4") * (-2)
     + (*ta1_)("x0, a4, c3, c5");
  madness::World::get_default().gop.fence();
}

void Task655::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c5, x1") * (*ta2_)("c5, c3, a2, x1");
  madness::World::get_default().gop.fence();
}

void Task656::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c3, a2, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a2");
  madness::World::get_default().gop.fence();
}

void Task657::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c5, c3, a2") += (*ta1_)("x0, c5, c3, a2") * (-2)
     + (*ta1_)("x0, a2, c3, c5");
  madness::World::get_default().gop.fence();
}

void Task658::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c5, x1") * (*ta2_)("c5, c3, a4, x1");
  madness::World::get_default().gop.fence();
}

void Task659::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c3, a4, x1") += (*ta1_)("x0, x1") * (*ta2_)("x0, c5, c3, a4");
  madness::World::get_default().gop.fence();
}

void Task660::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x0, c5, c3, a4") += (*ta1_)("x0, c5, c3, a4") * 4
     + (*ta1_)("x0, a4, c3, c5") * (-2);
  madness::World::get_default().gop.fence();
}

void Task661::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a4, a5, a2") * (*ta2_)("c1, a5, c3, x0");
  madness::World::get_default().gop.fence();
}

void Task662::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, a5, c3, x0") += (*ta1_)("x0, x1") * (*ta2_)("c1, a5, c3, x1") * 2;
  madness::World::get_default().gop.fence();
}

void Task663::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x0, a2, a5, a4") * (*ta2_)("c1, a5, c3, x0");
  madness::World::get_default().gop.fence();
}

void Task664::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, a5, c3, x0") += (*ta1_)("x0, x1") * (*ta2_)("c1, a5, c3, x1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task665::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a4, c1, x2") * (*ta2_)("c3, a2, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task666::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task667::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x1, x0") += (*ta1_)("c3, a2, x1, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a2") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task668::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("c3, x1, x0, a2") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task669::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x1, a2, c3, x0") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task670::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a2, c1, x2") * (*ta2_)("c3, a4, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task671::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a4, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task672::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x1, x0") += (*ta1_)("c3, a4, x1, x0") * (-1)
     + (*ta1_)("x1, a4, c3, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a4") * (-1);
  madness::World::get_default().gop.fence();
}

void Task673::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a4") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task674::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c3, c5") * (*ta2_)("a4, c5");
  madness::World::get_default().gop.fence();
}

void Task675::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, c5") += (*ta1_)("x1, x0") * (*ta2_)("x1, a4, c5, x0");
  madness::World::get_default().gop.fence();
}

void Task676::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, a4, c5, x0") += (*ta1_)("x1, a4, c5, x0") * 2
     + (*ta1_)("c5, a4, x1, x0") * (-4);
  madness::World::get_default().gop.fence();
}

void Task677::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c3, c5") * (*ta2_)("a2, c5");
  madness::World::get_default().gop.fence();
}

void Task678::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c5") += (*ta1_)("x1, x0") * (*ta2_)("x1, a2, c5, x0");
  madness::World::get_default().gop.fence();
}

void Task679::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, a2, c5, x0") += (*ta1_)("x1, a2, c5, x0") * (-1)
     + (*ta1_)("c5, a2, x1, x0") * 2;
  madness::World::get_default().gop.fence();
}

void Task680::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a4, a5, a2") * (*ta2_)("a5, c1");
  madness::World::get_default().gop.fence();
}

void Task681::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a5, c1, x0");
  madness::World::get_default().gop.fence();
}

void Task682::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, a5, c1, x0") += (*ta1_)("x1, a5, c1, x0") * (-2)
     + (*ta1_)("c1, a5, x1, x0") * 4;
  madness::World::get_default().gop.fence();
}

void Task683::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c3, a2, a5, a4") * (*ta2_)("a5, c1");
  madness::World::get_default().gop.fence();
}

void Task684::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c1") += (*ta1_)("x1, x0") * (*ta2_)("x1, a5, c1, x0");
  madness::World::get_default().gop.fence();
}

void Task685::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x1, a5, c1, x0") += (*ta1_)("x1, a5, c1, x0")
     + (*ta1_)("c1, a5, x1, x0") * (-2);
  madness::World::get_default().gop.fence();
}

void Task686::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, x3, x2") * (*ta2_)("c3, a2, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task687::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task688::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x1, x0") += (*ta1_)("c3, a2, x1, x0") * (-1)
     + (*ta1_)("x1, a2, c3, x0") * 0.5
     + (*ta1_)("x1, x0, c3, a2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task689::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a2, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a2") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task690::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, x3, x2") * (*ta2_)("c3, a4, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task691::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, a4, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task692::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x1, x0") += (*ta1_)("c3, a4, x1, x0") * 2
     + (*ta1_)("x1, a4, c3, x0") * (-1)
     + (*ta1_)("x1, x0, c3, a4") * 2;
  madness::World::get_default().gop.fence();
}

void Task693::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, a4, x3, x2") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, x1, x0, a4");
  madness::World::get_default().gop.fence();
}

void Task694::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("a4, x0");
  madness::World::get_default().gop.fence();
}

void Task695::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task696::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a4") * (*ta2_)("a2, x0");
  madness::World::get_default().gop.fence();
}

void Task697::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a2, x2, x1") * (-2);
  madness::World::get_default().gop.fence();
}

void Task698::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a4") * (*ta2_)("c3, c6, a5, a2") * (-8);
  madness::World::get_default().gop.fence();
}

void Task699::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c6, a5") * (*ta2_)("c3, c6, a5, a2") * 4;
  madness::World::get_default().gop.fence();
}

#endif
