//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks11.cc
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

#include <src/smith/RelMRCI_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task500::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, a5") += (*ta1_)("x1, x0") * (*ta2_)("c1, a5, x1, x0") * (-1);
  madness::World::get_default().gop.fence();
}

void Task501::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a2") * (*ta2_)("a4, x0");
  madness::World::get_default().gop.fence();
}

void Task502::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a4, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task503::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, x0, c3, a4") * (*ta2_)("a2, x0");
  madness::World::get_default().gop.fence();
}

void Task504::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a2, x2, x1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task505::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a4") * (*ta2_)("c3, c6, a5, a2") * (-2);
  madness::World::get_default().gop.fence();
}

void Task506::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c6, a5") * (*ta2_)("c3, c6, a5, a2") * 2;
  madness::World::get_default().gop.fence();
}

void Task507::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a5, c6, a2") * (*ta2_)("c3, c6, a5, a4") * 2;
  madness::World::get_default().gop.fence();
}

void Task508::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c6, a5") * (*ta2_)("c3, c6, a5, a4") * (-2);
  madness::World::get_default().gop.fence();
}

void Task509::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a6, c5, a4") * (*ta2_)("c3, a2, a6, c5") * 2;
  madness::World::get_default().gop.fence();
}

void Task510::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a4, c5, a6") * (*ta2_)("c3, a2, a6, c5") * (-2);
  madness::World::get_default().gop.fence();
}

void Task511::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a6, c5, a2") * (*ta2_)("c3, a4, a6, c5") * (-2);
  madness::World::get_default().gop.fence();
}

void Task512::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("c1, a2, c5, a6") * (*ta2_)("c3, a4, a6, c5") * 2;
  madness::World::get_default().gop.fence();
}

void Task513::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a4, c1, a2") * (*ta2_)("c3, x3");
  madness::World::get_default().gop.fence();
}

void Task514::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, x2, x1, x0") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task515::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x2, x1, c3, x0") * (-0.5);
  madness::World::get_default().gop.fence();
}

void Task516::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x3, a2, c1, a4") * (*ta2_)("c3, x3");
  madness::World::get_default().gop.fence();
}

void Task517::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c3, x2, x1, x0") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task518::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, x3") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x2, x1, c3, x0") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task519::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c5, a2") * (*ta2_)("c1, c3, c5, x1");
  madness::World::get_default().gop.fence();
}

void Task520::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c3, c5, x1") += (*ta1_)("x1, x0") * (*ta2_)("c1, x0, c3, c5") * (-1);
  madness::World::get_default().gop.fence();
}

void Task521::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c5, a4") * (*ta2_)("c1, c3, c5, x1");
  madness::World::get_default().gop.fence();
}

void Task522::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c3, c5, x1") += (*ta1_)("x1, x0") * (*ta2_)("c1, x0, c3, c5");
  madness::World::get_default().gop.fence();
}

void Task523::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a5, c1, a4") * (*ta2_)("a5, c3, a2, x1");
  madness::World::get_default().gop.fence();
}

void Task524::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c3, a2, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a2");
  madness::World::get_default().gop.fence();
}

void Task525::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, x0, c3, a2") += (*ta1_)("a5, x0, c3, a2") * (-1)
     + (*ta1_)("c3, x0, a5, a2");
  madness::World::get_default().gop.fence();
}

void Task526::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a4, c1, a5") * (*ta2_)("a5, c3, a2, x1");
  madness::World::get_default().gop.fence();
}

void Task527::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c3, a2, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a2");
  madness::World::get_default().gop.fence();
}

void Task528::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, x0, c3, a2") += (*ta1_)("a5, x0, c3, a2")
     + (*ta1_)("c3, x0, a5, a2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task529::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a5, c1, a2") * (*ta2_)("a5, c3, a4, x1");
  madness::World::get_default().gop.fence();
}

void Task530::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c3, a4, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a4");
  madness::World::get_default().gop.fence();
}

void Task531::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, x0, c3, a4") += (*ta1_)("a5, x0, c3, a4")
     + (*ta1_)("c3, x0, a5, a4") * (-1);
  madness::World::get_default().gop.fence();
}

void Task532::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a2, c1, a4, c3") += (*ta1_)("x1, a2, c1, a5") * (*ta2_)("a5, c3, a4, x1");
  madness::World::get_default().gop.fence();
}

void Task533::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, c3, a4, x1") += (*ta1_)("x1, x0") * (*ta2_)("a5, x0, c3, a4");
  madness::World::get_default().gop.fence();
}

void Task534::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("a5, x0, c3, a4") += (*ta1_)("a5, x0, c3, a4") * (-1)
     + (*ta1_)("c3, x0, a5, a4");
  madness::World::get_default().gop.fence();
}

void Task535::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, x0, a1");
  madness::World::get_default().gop.fence();
}

void Task536::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a1") * (*ta2_)("c2, a3, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task537::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a3, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task538::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("x1, a3") * (*ta2_)("c2, a1, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task539::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a1, x1, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("c2, a1, x3, x2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task540::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1") * (*ta2_)("a3, x0");
  madness::World::get_default().gop.fence();
}

void Task541::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a3, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a3, x2, x1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task542::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3") * (*ta2_)("a1, x0");
  madness::World::get_default().gop.fence();
}

void Task543::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("a1, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, a1, x2, x1");
  madness::World::get_default().gop.fence();
}

void Task544::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a3, c4, a1") * (*ta2_)("c4, x0");
  madness::World::get_default().gop.fence();
}

void Task545::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c4") * (-2);
  madness::World::get_default().gop.fence();
}

void Task546::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x0, x2, x1") * (*ta2_)("x3, c4, x2, x1") * (-1);
  madness::World::get_default().gop.fence();
}

void Task547::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c4, x0") += (*ta1_)("x3, x2, x1, x0") * (*ta2_)("x3, x2, x1, c4") * (-1);
  madness::World::get_default().gop.fence();
}

void Task548::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, a3, x0, a1") += (*ta1_)("c2, a1, c4, a3") * (*ta2_)("c4, x0");
  madness::World::get_default().gop.fence();
}

void Task549::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c4, x0") += (*ta1_)("x1, x0") * (*ta2_)("x1, c4") * 2;
  madness::World::get_default().gop.fence();
}

#endif
