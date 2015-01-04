//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_deciqq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_deciq() {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor62 = {deci};
  auto task62 = make_shared<Task62>(tensor62);
  deciq->add_task(task62);

  vector<IndexRange> I62_index = {ci_};
  auto I62 = make_shared<Tensor>(I62_index);
  vector<shared_ptr<Tensor>> tensor63 = {deci, I62};
  auto task63 = make_shared<Task63>(tensor63, cindex);
  task63->add_dep(task62);
  deciq->add_task(task63);

  vector<IndexRange> I63_index;
  auto I63 = make_shared<Tensor>(I63_index);
  vector<shared_ptr<Tensor>> tensor64 = {I62, Gamma6_(), I63};
  auto task64 = make_shared<Task64>(tensor64, cindex);
  task63->add_dep(task64);
  task64->add_dep(task62);
  deciq->add_task(task64);

  vector<IndexRange> I64_index = {virt_, closed_, virt_, closed_};
  auto I64 = make_shared<Tensor>(I64_index);
  vector<shared_ptr<Tensor>> tensor65 = {I63, t2, I64};
  auto task65 = make_shared<Task65>(tensor65, cindex);
  task64->add_dep(task65);
  task65->add_dep(task62);
  deciq->add_task(task65);

  vector<shared_ptr<Tensor>> tensor66 = {I64, t2};
  auto task66 = make_shared<Task66>(tensor66, cindex);
  task65->add_dep(task66);
  task66->add_dep(task62);
  deciq->add_task(task66);

  vector<IndexRange> I67_index = {closed_, virt_, closed_, virt_};
  auto I67 = make_shared<Tensor>(I67_index);
  vector<shared_ptr<Tensor>> tensor67 = {I63, t2, I67};
  auto task67 = make_shared<Task67>(tensor67, cindex);
  task64->add_dep(task67);
  task67->add_dep(task62);
  deciq->add_task(task67);

  vector<shared_ptr<Tensor>> tensor68 = {I67, t2};
  auto task68 = make_shared<Task68>(tensor68, cindex);
  task67->add_dep(task68);
  task68->add_dep(task62);
  deciq->add_task(task68);

  vector<IndexRange> I73_index = {virt_, closed_, virt_, closed_};
  auto I73 = make_shared<Tensor>(I73_index);
  vector<shared_ptr<Tensor>> tensor69 = {I63, t2, I73};
  auto task69 = make_shared<Task69>(tensor69, cindex);
  task64->add_dep(task69);
  task69->add_dep(task62);
  deciq->add_task(task69);

  vector<shared_ptr<Tensor>> tensor70 = {I73, t2};
  auto task70 = make_shared<Task70>(tensor70, cindex);
  task69->add_dep(task70);
  task70->add_dep(task62);
  deciq->add_task(task70);

  return deciq;
}


