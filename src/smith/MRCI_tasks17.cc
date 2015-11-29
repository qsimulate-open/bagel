//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks17.cc
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

#include <src/smith/MRCI_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task800::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, c4, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a1, c4, x2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task801::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, c4, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c4, a1, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task802::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a3, a4, a1") * (*ta2_)("a4, c2, x0, x1");
  madness::World::get_default().gop.fence();
}

void Task803::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, c2, x0, x1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("x3, a4, c2, x2");
  madness::World::get_default().gop.fence();
}

void Task804::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, c2, x0, x1") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a4, x3, x2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task805::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x1, a1, a4, a3") * (*ta2_)("a4, c2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task806::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, c2, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, a4, c2, x2");
  madness::World::get_default().gop.fence();
}

void Task807::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, a4, c2, x2") += (*ta1_)("x3, a4, c2, x2") * (-1)
     + (*ta1_)("c2, a4, x3, x2") * 2;
  madness::World::get_default().gop.fence();
}

void Task808::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a1, x5, x4, x0");
  madness::World::get_default().gop.fence();
}

void Task809::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task810::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, a1");
  madness::World::get_default().gop.fence();
}

void Task811::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a1, x5, x4") * (*ta2_)("a3, x5, x4, x0");
  madness::World::get_default().gop.fence();
}

void Task812::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x5, x4, x0") += (*ta1_)("x5, x4, x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task813::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x5, x4, x0") += (*ta1_)("x5, x4, x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, a3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task814::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a1, x2, x1") * (*ta2_)("a3, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task815::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task816::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a3, x2, x1") * (*ta2_)("a1, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task817::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task818::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, x2, x1, a1") * (*ta2_)("a3, x2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task819::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x2, x1, x0") += (*ta1_)("x5, x2, x4, x3, x1, x0") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task820::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, x2, x1, a3") * (*ta2_)("a1, x0, x1, x2");
  madness::World::get_default().gop.fence();
}

void Task821::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x1, x2") += (*ta1_)("x5, x0, x4, x3, x1, x2") * (*ta2_)("x5, a1, x4, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task822::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x2, a1, c2, x1") * (*ta2_)("a3, x1, x2, x0");
  madness::World::get_default().gop.fence();
}

void Task823::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x1, x2, x0") += (*ta1_)("x5, x1, x4, x3, x2, x0") * (*ta2_)("x5, a3, x4, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task824::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x2, a3, c2, x1") * (*ta2_)("a1, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task825::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task826::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x2, x1, c2, a1") * (*ta2_)("a3, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task827::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a3, x4, x3") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task828::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x2, x1, c2, a3") * (*ta2_)("a1, x0, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task829::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0, x2, x1") += (*ta1_)("x5, x0, x4, x3, x2, x1") * (*ta2_)("x5, a1, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task830::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a3, a4, a1") * (*ta2_)("a4, x0");
  madness::World::get_default().gop.fence();
}

void Task831::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1") * 2;
  madness::World::get_default().gop.fence();
}

void Task832::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c2, a1, a4, a3") * (*ta2_)("a4, x0");
  madness::World::get_default().gop.fence();
}

void Task833::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task834::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c4, a3, c5, a1") * (*ta2_)("c5, c2, c4, x0");
  madness::World::get_default().gop.fence();
}

void Task835::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c2, c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c5, c2, c4") * 4;
  madness::World::get_default().gop.fence();
}

void Task836::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("c4, a1, c5, a3") * (*ta2_)("c5, c2, c4, x0");
  madness::World::get_default().gop.fence();
}

void Task837::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c5, c2, c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c5, c2, c4") * (-2);
  madness::World::get_default().gop.fence();
}

void Task838::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, a3, c4, a1") * (*ta2_)("c2, c4, x3, x0");
  madness::World::get_default().gop.fence();
}

void Task839::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("c2, c4, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task840::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x2, x1") += (*ta1_)("c2, c4, x2, x1") * 0.5
     + (*ta1_)("x2, x1, c2, c4") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task841::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, x2, x1, c4") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task842::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x1, x2, x0") * (*ta2_)("x2, c4, c2, x1") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task843::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, a1, c4, a3") * (*ta2_)("c2, c4, x3, x0");
  madness::World::get_default().gop.fence();
}

void Task844::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("c2, c4, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task845::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x2, x1") += (*ta1_)("c2, c4, x2, x1") * (-1)
     + (*ta1_)("x2, c4, c2, x1") * 0.5
     + (*ta1_)("x2, x1, c2, c4") * (-1);
  madness::World::get_default().gop.fence();
}

void Task846::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c4, x3, x0") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("c2, x2, x1, c4") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task847::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, c2, x0, a1") += (*ta1_)("x3, a4, c2, a3") * (*ta2_)("a4, a1, x3, x0");
  madness::World::get_default().gop.fence();
}

void Task848::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, a1, x3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("a4, a1, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task849::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, a1, x2, x1") += (*ta1_)("a4, a1, x2, x1")
     + (*ta1_)("x2, x1, a4, a1");
  madness::World::get_default().gop.fence();
}

#endif
