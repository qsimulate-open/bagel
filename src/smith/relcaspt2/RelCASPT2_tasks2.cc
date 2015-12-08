//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_tasks2.cc
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

#include <src/smith/relcaspt2/RelCASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASPT2;

void Task50::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("a2, x3");
}

void Task51::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
}

void Task52::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, x3");
}

void Task53::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("a3, c2");
}

void Task54::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, c2") * (-1);
}

void Task55::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x3, x2, x4, x1, x0") * (*ta2_)("x3, x5, c1, x4");
}

void Task56::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("a2, x3") * (-1);
}

void Task57::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3") * 0.5;
}

void Task58::compute_() {
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c1, c3, x0, a2");
}

void Task59::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x1, a2") * (*ta2_)("c1, c3, x1, x0");
}

void Task60::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
}

void Task61::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2") * (*ta2_)("c3, x0");
}

void Task62::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c3, x1");
}

void Task63::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, a2") * (*ta2_)("c1, x0");
}

void Task64::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, x2, c1, x1") * (-1);
}

void Task65::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x0, x3") * (*ta2_)("c3, a2, c1, x3");
}

void Task66::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c3, a2, c1, x3") += (*ta1_)("c3, a2, c1, x3") * (-1)
     + (*ta1_)("c1, a2, c3, x3");
}

void Task67::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x0, x1") * (*ta2_)("c1, a2, c3, x1");
}

void Task68::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c3, a2, c1, x1") * e0_
     + (*ta1_)("c1, a2, c3, x1") * e0_ * (-1)
     + (*ta2_)("c3, x1, c1, a2")
     + (*ta2_)("c1, x1, c3, a2") * (-1);
}

void Task69::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c4, a2, c3, x1") * (*ta2_)("c1, c4") * (-1);
}

void Task70::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c3, a2, c4, x1") * (*ta2_)("c1, c4");
}

void Task71::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c4, a2, c1, x1") * (*ta2_)("c3, c4");
}

void Task72::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c3, a4, c1, x1") * (*ta2_)("a4, a2") * (-1);
}

void Task73::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c1, a2, c4, x1") * (*ta2_)("c3, c4") * (-1);
}

void Task74::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, a2, c3, x1") += (*ta1_)("c1, a4, c3, x1") * (*ta2_)("a4, a2");
}

void Task75::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, x1") * (*ta2_)("a2, c3, x1, x0");
}

void Task76::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x1, x0, x2") * (*ta2_)("x3, a2, c3, x2");
}

void Task77::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c3, x1, x0") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("c3, a2, x3, x2") * (-1);
}

void Task78::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c3, x1") * (*ta2_)("a2, c1, x0, x1");
}

void Task79::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a2, c1, x0, x1") += (*ta1_)("x3, x2, x0, x1") * (*ta2_)("x3, a2, c1, x2");
}

void Task80::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, a2, c1, x2") += (*ta1_)("x3, a2, c1, x2") * (-1)
     + (*ta1_)("c1, a2, x3, x2");
}

void Task81::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a4, c3, a2") * (*ta2_)("a4, x0");
}

void Task82::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x0, x1") * (*ta2_)("a4, x1") * (-2);
}

void Task83::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2, c3, a4") * (*ta2_)("a4, x0");
}

void Task84::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a4, x0") += (*ta1_)("x0, x1") * (*ta2_)("a4, x1") * 2;
}

void Task85::compute_() {
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("c2, x1, x0, a1");
}

void Task86::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x2, a1") * (*ta2_)("c2, x1, x2, x0");
}

void Task87::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x2, x0") += (*ta1_)("x5, x4, x1, x3, x2, x0") * (*ta2_)("x5, x4, c2, x3") * (-1);
}

void Task88::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x1, x3, x2, x0") * (*ta2_)("x2, a1, c2, x3");
}

void Task89::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1, c2, x3") += (*ta1_)("c2, x3, x2, a1") * (-0.5);
}

void Task90::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, a1, c2, x3") += (*ta1_)("c3, a1, c2, x3") * (*ta2_)("x2, c3");
}

void Task91::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("x2, c2, a1, x3");
}

void Task92::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c2, a1, x3") += (*ta1_)("c2, a1, c3, x3") * (*ta2_)("x2, c3") * (-1);
}

void Task93::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x0, x1, x4") * (*ta2_)("x5, a1, c2, x4");
}

void Task94::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x3, x0, x1, x2") * (*ta2_)("c2, x3, a1, x2");
}

void Task95::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, a1, x2") += (*ta1_)("x3, a1, c2, x2") * e0_ * (-1)
     + (*ta2_)("x3, a1, c2, x2") * 0.5;
}

void Task96::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, a1, x2") += (*ta1_)("x3, a1, c3, x2") * (*ta2_)("c2, c3") * (-1);
}

void Task97::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, a1, x2") += (*ta1_)("x3, a3, c2, x2") * (*ta2_)("a3, a1");
}

void Task98::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, a1, x2") += (*ta1_)("x3, a1, c2, a3") * (*ta2_)("a3, x2");
}

void Task99::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x1, x0, a1") += (*ta1_)("x5, x4, x1, x0") * (*ta2_)("c2, a1, x5, x4") * (-1);
}

#endif
