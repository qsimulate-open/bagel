//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density2qq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_density2q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor1065 = {Den1};
  auto task1065 = make_shared<Task1065>(tensor1065);
  density2q->add_task(task1065);

  vector<IndexRange> I1170_index = {closed_, closed_, active_, active_};
  auto I1170 = make_shared<Tensor>(I1170_index);
  vector<shared_ptr<Tensor>> tensor1066 = {Den1, I1170};
  auto task1066 = make_shared<Task1066>(tensor1066, pindex);
  task1066->add_dep(task1065);
  density2q->add_task(task1066);

  vector<IndexRange> I1171_index = {closed_, active_, closed_, active_};
  auto I1171 = make_shared<Tensor>(I1171_index);
  vector<shared_ptr<Tensor>> tensor1067 = {I1170, Gamma94_(), I1171};
  auto task1067 = make_shared<Task1067>(tensor1067, pindex);
  task1066->add_dep(task1067);
  task1067->add_dep(task1065);
  density2q->add_task(task1067);

  vector<shared_ptr<Tensor>> tensor1068 = {I1171, t2};
  auto task1068 = make_shared<Task1068>(tensor1068, pindex);
  task1067->add_dep(task1068);
  task1068->add_dep(task1065);
  density2q->add_task(task1068);

  vector<IndexRange> I1172_index = {closed_, active_, active_, active_};
  auto I1172 = make_shared<Tensor>(I1172_index);
  vector<shared_ptr<Tensor>> tensor1069 = {Den1, I1172};
  auto task1069 = make_shared<Task1069>(tensor1069, pindex);
  task1069->add_dep(task1065);
  density2q->add_task(task1069);

  vector<IndexRange> I1173_index = {active_, active_, closed_, active_};
  auto I1173 = make_shared<Tensor>(I1173_index);
  vector<shared_ptr<Tensor>> tensor1070 = {I1172, Gamma6_(), I1173};
  auto task1070 = make_shared<Task1070>(tensor1070, pindex);
  task1069->add_dep(task1070);
  task1070->add_dep(task1065);
  density2q->add_task(task1070);

  vector<shared_ptr<Tensor>> tensor1071 = {I1173, t2};
  auto task1071 = make_shared<Task1071>(tensor1071, pindex);
  task1070->add_dep(task1071);
  task1071->add_dep(task1065);
  density2q->add_task(task1071);

  vector<IndexRange> I1174_index = {closed_, virt_, closed_, active_};
  auto I1174 = make_shared<Tensor>(I1174_index);
  vector<shared_ptr<Tensor>> tensor1072 = {Den1, I1174};
  auto task1072 = make_shared<Task1072>(tensor1072, pindex);
  task1072->add_dep(task1065);
  density2q->add_task(task1072);

  vector<IndexRange> I1175_index = {closed_, virt_, closed_, active_};
  auto I1175 = make_shared<Tensor>(I1175_index);
  vector<shared_ptr<Tensor>> tensor1073 = {I1174, Gamma16_(), I1175};
  auto task1073 = make_shared<Task1073>(tensor1073, pindex);
  task1072->add_dep(task1073);
  task1073->add_dep(task1065);
  density2q->add_task(task1073);

  vector<shared_ptr<Tensor>> tensor1074 = {I1175, t2};
  auto task1074 = make_shared<Task1074>(tensor1074, pindex);
  task1073->add_dep(task1074);
  task1074->add_dep(task1065);
  density2q->add_task(task1074);

  vector<IndexRange> I1178_index = {virt_, closed_, active_, active_};
  auto I1178 = make_shared<Tensor>(I1178_index);
  vector<shared_ptr<Tensor>> tensor1075 = {Den1, I1178};
  auto task1075 = make_shared<Task1075>(tensor1075, pindex);
  task1075->add_dep(task1065);
  density2q->add_task(task1075);

  vector<IndexRange> I1179_index = {active_, virt_, closed_, active_};
  auto I1179 = make_shared<Tensor>(I1179_index);
  vector<shared_ptr<Tensor>> tensor1076 = {I1178, Gamma32_(), I1179};
  auto task1076 = make_shared<Task1076>(tensor1076, pindex);
  task1075->add_dep(task1076);
  task1076->add_dep(task1065);
  density2q->add_task(task1076);

  vector<shared_ptr<Tensor>> tensor1077 = {I1179, t2};
  auto task1077 = make_shared<Task1077>(tensor1077, pindex);
  task1076->add_dep(task1077);
  task1077->add_dep(task1065);
  density2q->add_task(task1077);

  vector<IndexRange> I1181_index = {closed_, virt_, active_, active_};
  auto I1181 = make_shared<Tensor>(I1181_index);
  vector<shared_ptr<Tensor>> tensor1078 = {I1178, Gamma35_(), I1181};
  auto task1078 = make_shared<Task1078>(tensor1078, pindex);
  task1075->add_dep(task1078);
  task1078->add_dep(task1065);
  density2q->add_task(task1078);

  vector<shared_ptr<Tensor>> tensor1079 = {I1181, t2};
  auto task1079 = make_shared<Task1079>(tensor1079, pindex);
  task1078->add_dep(task1079);
  task1079->add_dep(task1065);
  density2q->add_task(task1079);

  vector<IndexRange> I1182_index = {virt_, closed_, active_, active_};
  auto I1182 = make_shared<Tensor>(I1182_index);
  vector<shared_ptr<Tensor>> tensor1080 = {Den1, I1182};
  auto task1080 = make_shared<Task1080>(tensor1080, pindex);
  task1080->add_dep(task1065);
  density2q->add_task(task1080);

  vector<IndexRange> I1183_index = {active_, virt_, closed_, active_};
  auto I1183 = make_shared<Tensor>(I1183_index);
  vector<shared_ptr<Tensor>> tensor1081 = {I1182, Gamma35_(), I1183};
  auto task1081 = make_shared<Task1081>(tensor1081, pindex);
  task1080->add_dep(task1081);
  task1081->add_dep(task1065);
  density2q->add_task(task1081);

  vector<shared_ptr<Tensor>> tensor1082 = {I1183, t2};
  auto task1082 = make_shared<Task1082>(tensor1082, pindex);
  task1081->add_dep(task1082);
  task1082->add_dep(task1065);
  density2q->add_task(task1082);

  vector<IndexRange> I1186_index = {virt_, active_, active_, active_};
  auto I1186 = make_shared<Tensor>(I1186_index);
  vector<shared_ptr<Tensor>> tensor1083 = {Den1, I1186};
  auto task1083 = make_shared<Task1083>(tensor1083, pindex);
  task1083->add_dep(task1065);
  density2q->add_task(task1083);

  vector<IndexRange> I1187_index = {active_, virt_, active_, active_};
  auto I1187 = make_shared<Tensor>(I1187_index);
  vector<shared_ptr<Tensor>> tensor1084 = {I1186, Gamma59_(), I1187};
  auto task1084 = make_shared<Task1084>(tensor1084, pindex);
  task1083->add_dep(task1084);
  task1084->add_dep(task1065);
  density2q->add_task(task1084);

  vector<shared_ptr<Tensor>> tensor1085 = {I1187, t2};
  auto task1085 = make_shared<Task1085>(tensor1085, pindex);
  task1084->add_dep(task1085);
  task1085->add_dep(task1065);
  density2q->add_task(task1085);

  vector<IndexRange> I1188_index = {closed_, virt_, closed_, virt_};
  auto I1188 = make_shared<Tensor>(I1188_index);
  vector<shared_ptr<Tensor>> tensor1086 = {Den1, I1188};
  auto task1086 = make_shared<Task1086>(tensor1086, pindex);
  task1086->add_dep(task1065);
  density2q->add_task(task1086);

  vector<shared_ptr<Tensor>> tensor1087 = {I1188, t2};
  auto task1087 = make_shared<Task1087>(tensor1087, pindex);
  task1086->add_dep(task1087);
  task1087->add_dep(task1065);
  density2q->add_task(task1087);

  vector<IndexRange> I1190_index = {virt_, closed_, virt_, active_};
  auto I1190 = make_shared<Tensor>(I1190_index);
  vector<shared_ptr<Tensor>> tensor1088 = {Den1, I1190};
  auto task1088 = make_shared<Task1088>(tensor1088, pindex);
  task1088->add_dep(task1065);
  density2q->add_task(task1088);

  vector<IndexRange> I1191_index = {active_, virt_, closed_, virt_};
  auto I1191 = make_shared<Tensor>(I1191_index);
  vector<shared_ptr<Tensor>> tensor1089 = {I1190, Gamma38_(), I1191};
  auto task1089 = make_shared<Task1089>(tensor1089, pindex);
  task1088->add_dep(task1089);
  task1089->add_dep(task1065);
  density2q->add_task(task1089);

  vector<shared_ptr<Tensor>> tensor1090 = {I1191, t2};
  auto task1090 = make_shared<Task1090>(tensor1090, pindex);
  task1089->add_dep(task1090);
  task1090->add_dep(task1065);
  density2q->add_task(task1090);

  vector<IndexRange> I1194_index = {virt_, virt_, active_, active_};
  auto I1194 = make_shared<Tensor>(I1194_index);
  vector<shared_ptr<Tensor>> tensor1091 = {Den1, I1194};
  auto task1091 = make_shared<Task1091>(tensor1091, pindex);
  task1091->add_dep(task1065);
  density2q->add_task(task1091);

  vector<IndexRange> I1195_index = {active_, virt_, active_, virt_};
  auto I1195 = make_shared<Tensor>(I1195_index);
  vector<shared_ptr<Tensor>> tensor1092 = {I1194, Gamma60_(), I1195};
  auto task1092 = make_shared<Task1092>(tensor1092, pindex);
  task1091->add_dep(task1092);
  task1092->add_dep(task1065);
  density2q->add_task(task1092);

  vector<shared_ptr<Tensor>> tensor1093 = {I1195, t2};
  auto task1093 = make_shared<Task1093>(tensor1093, pindex);
  task1092->add_dep(task1093);
  task1093->add_dep(task1065);
  density2q->add_task(task1093);

  return density2q;
}


