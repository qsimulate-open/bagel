//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks18.h>
#include <src/smith/mrci/MRCI_tasks19.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq9(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I239_index = {virt_, active_, active_, virt_};
  auto I239 = make_shared<Tensor>(I239_index);
  auto tensor880 = vector<shared_ptr<Tensor>>{r, I239};
  auto task880 = make_shared<Task880>(tensor880, pindex);
  task880->add_dep(task108);
  residualq->add_task(task880);

  vector<IndexRange> I240_index = {virt_, active_, active_, active_};
  auto I240 = make_shared<Tensor>(I240_index);
  auto tensor881 = vector<shared_ptr<Tensor>>{I239, h1_, I240};
  auto task881 = make_shared<Task881>(tensor881, pindex);
  task880->add_dep(task881);
  task881->add_dep(task108);
  residualq->add_task(task881);

  auto tensor882 = vector<shared_ptr<Tensor>>{I240, Gamma50_(), t2};
  auto task882 = make_shared<Task882>(tensor882, pindex);
  task881->add_dep(task882);
  task882->add_dep(task108);
  residualq->add_task(task882);

  vector<IndexRange> I243_index = {active_, active_, virt_, virt_};
  auto I243 = make_shared<Tensor>(I243_index);
  auto tensor883 = vector<shared_ptr<Tensor>>{I239, Gamma51_(), I243};
  auto task883 = make_shared<Task883>(tensor883, pindex);
  task880->add_dep(task883);
  task883->add_dep(task108);
  residualq->add_task(task883);

  auto tensor884 = vector<shared_ptr<Tensor>>{I243, t2, h1_};
  auto task884 = make_shared<Task884>(tensor884, pindex);
  task883->add_dep(task884);
  task884->add_dep(task108);
  residualq->add_task(task884);

  auto tensor885 = vector<shared_ptr<Tensor>>{I243, t2, h1_};
  auto task885 = make_shared<Task885>(tensor885, pindex);
  task883->add_dep(task885);
  task885->add_dep(task108);
  residualq->add_task(task885);

  auto tensor886 = vector<shared_ptr<Tensor>>{I243, t2, v2_};
  auto task886 = make_shared<Task886>(tensor886, pindex);
  task883->add_dep(task886);
  task886->add_dep(task108);
  residualq->add_task(task886);

  auto tensor887 = vector<shared_ptr<Tensor>>{I243, t2, v2_};
  auto task887 = make_shared<Task887>(tensor887, pindex);
  task883->add_dep(task887);
  task887->add_dep(task108);
  residualq->add_task(task887);

  auto tensor888 = vector<shared_ptr<Tensor>>{I243, t2, v2_};
  auto task888 = make_shared<Task888>(tensor888, pindex);
  task883->add_dep(task888);
  task888->add_dep(task108);
  residualq->add_task(task888);

  vector<IndexRange> I1622_index = {virt_, closed_, active_, active_, active_, active_};
  auto I1622 = make_shared<Tensor>(I1622_index);
  auto tensor889 = vector<shared_ptr<Tensor>>{I239, t2, I1622};
  auto task889 = make_shared<Task889>(tensor889, pindex);
  task880->add_dep(task889);
  task889->add_dep(task108);
  residualq->add_task(task889);

  auto tensor890 = vector<shared_ptr<Tensor>>{I1622, Gamma531_(), v2_};
  auto task890 = make_shared<Task890>(tensor890, pindex);
  task889->add_dep(task890);
  task890->add_dep(task108);
  residualq->add_task(task890);

  vector<IndexRange> I1625_index = {virt_, closed_, active_, active_, active_, active_};
  auto I1625 = make_shared<Tensor>(I1625_index);
  auto tensor891 = vector<shared_ptr<Tensor>>{I239, t2, I1625};
  auto task891 = make_shared<Task891>(tensor891, pindex);
  task880->add_dep(task891);
  task891->add_dep(task108);
  residualq->add_task(task891);

  auto tensor892 = vector<shared_ptr<Tensor>>{I1625, Gamma532_(), v2_};
  auto task892 = make_shared<Task892>(tensor892, pindex);
  task891->add_dep(task892);
  task892->add_dep(task108);
  residualq->add_task(task892);

  vector<IndexRange> I1628_index = {virt_, active_, active_, active_, active_, active_};
  auto I1628 = make_shared<Tensor>(I1628_index);
  auto tensor893 = vector<shared_ptr<Tensor>>{I239, t2, I1628};
  auto task893 = make_shared<Task893>(tensor893, pindex);
  task880->add_dep(task893);
  task893->add_dep(task108);
  residualq->add_task(task893);

  auto tensor894 = vector<shared_ptr<Tensor>>{I1628, Gamma533_(), v2_};
  auto task894 = make_shared<Task894>(tensor894, pindex);
  task893->add_dep(task894);
  task894->add_dep(task108);
  residualq->add_task(task894);

  auto tensor895 = vector<shared_ptr<Tensor>>{I1628, Gamma349_(), v2_};
  auto task895 = make_shared<Task895>(tensor895, pindex);
  task893->add_dep(task895);
  task895->add_dep(task108);
  residualq->add_task(task895);

  vector<IndexRange> I1634_index = {virt_, active_, active_, active_};
  auto I1634 = make_shared<Tensor>(I1634_index);
  auto tensor896 = vector<shared_ptr<Tensor>>{I239, v2_, I1634};
  auto task896 = make_shared<Task896>(tensor896, pindex);
  task880->add_dep(task896);
  task896->add_dep(task108);
  residualq->add_task(task896);

  auto tensor897 = vector<shared_ptr<Tensor>>{I1634, Gamma471_(), t2};
  auto task897 = make_shared<Task897>(tensor897, pindex);
  task896->add_dep(task897);
  task897->add_dep(task108);
  residualq->add_task(task897);

  vector<IndexRange> I1640_index = {closed_, active_, active_, active_};
  auto I1640 = make_shared<Tensor>(I1640_index);
  auto tensor898 = vector<shared_ptr<Tensor>>{I239, t2, I1640};
  auto task898 = make_shared<Task898>(tensor898, pindex);
  task880->add_dep(task898);
  task898->add_dep(task108);
  residualq->add_task(task898);

  auto tensor899 = vector<shared_ptr<Tensor>>{I1640, Gamma526_(), v2_};
  auto task899 = make_shared<Task899>(tensor899, pindex);
  task898->add_dep(task899);
  task899->add_dep(task108);
  residualq->add_task(task899);

  auto tensor900 = vector<shared_ptr<Tensor>>{I1640, Gamma50_(), v2_};
  auto task900 = make_shared<Task900>(tensor900, pindex);
  task898->add_dep(task900);
  task900->add_dep(task108);
  residualq->add_task(task900);

  vector<IndexRange> I1646_index = {active_, virt_, active_, virt_};
  auto I1646 = make_shared<Tensor>(I1646_index);
  auto tensor901 = vector<shared_ptr<Tensor>>{I239, Gamma503_(), I1646};
  auto task901 = make_shared<Task901>(tensor901, pindex);
  task880->add_dep(task901);
  task901->add_dep(task108);
  residualq->add_task(task901);

  auto tensor902 = vector<shared_ptr<Tensor>>{I1646, t2, v2_};
  auto task902 = make_shared<Task902>(tensor902, pindex);
  task901->add_dep(task902);
  task902->add_dep(task108);
  residualq->add_task(task902);

  vector<IndexRange> I1658_index = {virt_, active_, active_, active_, virt_, active_};
  auto I1658 = make_shared<Tensor>(I1658_index);
  auto tensor903 = vector<shared_ptr<Tensor>>{I239, Gamma526_(), I1658};
  auto task903 = make_shared<Task903>(tensor903, pindex);
  task880->add_dep(task903);
  task903->add_dep(task108);
  residualq->add_task(task903);

  vector<IndexRange> I1659_index = {virt_, virt_, active_, active_};
  auto I1659 = make_shared<Tensor>(I1659_index);
  auto tensor904 = vector<shared_ptr<Tensor>>{I1658, t2, I1659};
  auto task904 = make_shared<Task904>(tensor904, pindex);
  task903->add_dep(task904);
  task904->add_dep(task108);
  residualq->add_task(task904);

  auto tensor905 = vector<shared_ptr<Tensor>>{I1659, v2_};
  auto task905 = make_shared<Task905>(tensor905, pindex);
  task904->add_dep(task905);
  task905->add_dep(task108);
  residualq->add_task(task905);

  vector<IndexRange> I1661_index = {active_, active_, virt_, active_, virt_, active_};
  auto I1661 = make_shared<Tensor>(I1661_index);
  auto tensor906 = vector<shared_ptr<Tensor>>{I239, Gamma50_(), I1661};
  auto task906 = make_shared<Task906>(tensor906, pindex);
  task880->add_dep(task906);
  task906->add_dep(task108);
  residualq->add_task(task906);

  auto tensor907 = vector<shared_ptr<Tensor>>{I1661, t2, v2_};
  auto task907 = make_shared<Task907>(tensor907, pindex);
  task906->add_dep(task907);
  task907->add_dep(task108);
  residualq->add_task(task907);

  vector<IndexRange> I1664_index = {active_, virt_, active_, active_, virt_, active_};
  auto I1664 = make_shared<Tensor>(I1664_index);
  auto tensor908 = vector<shared_ptr<Tensor>>{I239, Gamma545_(), I1664};
  auto task908 = make_shared<Task908>(tensor908, pindex);
  task880->add_dep(task908);
  task908->add_dep(task108);
  residualq->add_task(task908);

  auto tensor909 = vector<shared_ptr<Tensor>>{I1664, t2, v2_};
  auto task909 = make_shared<Task909>(tensor909, pindex);
  task908->add_dep(task909);
  task909->add_dep(task108);
  residualq->add_task(task909);

  vector<IndexRange> I260_index = {closed_, closed_, active_, active_};
  auto I260 = make_shared<Tensor>(I260_index);
  auto tensor910 = vector<shared_ptr<Tensor>>{r, I260};
  auto task910 = make_shared<Task910>(tensor910, pindex);
  task910->add_dep(task108);
  residualq->add_task(task910);

  vector<IndexRange> I261_index = {closed_, closed_, active_, active_};
  auto I261 = make_shared<Tensor>(I261_index);
  auto tensor911 = vector<shared_ptr<Tensor>>{I260, Gamma2_(), I261};
  auto task911 = make_shared<Task911>(tensor911, pindex);
  task910->add_dep(task911);
  task911->add_dep(task108);
  residualq->add_task(task911);

  auto tensor912 = vector<shared_ptr<Tensor>>{I261, t2, v2_};
  auto task912 = make_shared<Task912>(tensor912, pindex);
  task911->add_dep(task912);
  task912->add_dep(task108);
  residualq->add_task(task912);

  auto tensor913 = vector<shared_ptr<Tensor>>{I261, t2, v2_};
  auto task913 = make_shared<Task913>(tensor913, pindex);
  task911->add_dep(task913);
  task913->add_dep(task108);
  residualq->add_task(task913);

  auto tensor914 = vector<shared_ptr<Tensor>>{I260, Gamma548_(), t2};
  auto task914 = make_shared<Task914>(tensor914, pindex);
  task910->add_dep(task914);
  task914->add_dep(task108);
  residualq->add_task(task914);

  auto tensor915 = vector<shared_ptr<Tensor>>{I260, Gamma549_(), t2};
  auto task915 = make_shared<Task915>(tensor915, pindex);
  task910->add_dep(task915);
  task915->add_dep(task108);
  residualq->add_task(task915);

  vector<IndexRange> I1112_index = {closed_, closed_, virt_, virt_};
  auto I1112 = make_shared<Tensor>(I1112_index);
  auto tensor916 = vector<shared_ptr<Tensor>>{r, I1112};
  auto task916 = make_shared<Task916>(tensor916, pindex);
  task916->add_dep(task108);
  residualq->add_task(task916);

  vector<IndexRange> I1113_index = {closed_, closed_, active_, active_};
  auto I1113 = make_shared<Tensor>(I1113_index);
  auto tensor917 = vector<shared_ptr<Tensor>>{I1112, v2_, I1113};
  auto task917 = make_shared<Task917>(tensor917, pindex);
  task916->add_dep(task917);
  task917->add_dep(task108);
  residualq->add_task(task917);

  auto tensor918 = vector<shared_ptr<Tensor>>{I1113, Gamma2_(), t2};
  auto task918 = make_shared<Task918>(tensor918, pindex);
  task917->add_dep(task918);
  task918->add_dep(task108);
  residualq->add_task(task918);

  shared_ptr<Task919> task919;
  if (diagonal) {
    auto tensor919 = vector<shared_ptr<Tensor>>{I1112, t2, v2_};
    task919 = make_shared<Task919>(tensor919, pindex);
    task916->add_dep(task919);
    task919->add_dep(task108);
    residualq->add_task(task919);
  }

  shared_ptr<Task920> task920;
  if (diagonal) {
    auto tensor920 = vector<shared_ptr<Tensor>>{I1112, t2, v2_};
    task920 = make_shared<Task920>(tensor920, pindex);
    task916->add_dep(task920);
    task920->add_dep(task108);
    residualq->add_task(task920);
  }

  vector<IndexRange> I1352_index = {closed_, closed_, active_, active_};
  auto I1352 = make_shared<Tensor>(I1352_index);
  auto tensor921 = vector<shared_ptr<Tensor>>{I1112, t2, I1352};
  auto task921 = make_shared<Task921>(tensor921, pindex);
  task916->add_dep(task921);
  task921->add_dep(task108);
  residualq->add_task(task921);

  auto tensor922 = vector<shared_ptr<Tensor>>{I1352, Gamma503_(), v2_};
  auto task922 = make_shared<Task922>(tensor922, pindex);
  task921->add_dep(task922);
  task922->add_dep(task108);
  residualq->add_task(task922);

  vector<IndexRange> I1693_index = {closed_, virt_, closed_, virt_};
  auto I1693 = make_shared<Tensor>(I1693_index);
  auto tensor923 = vector<shared_ptr<Tensor>>{I1112, Gamma558_(), I1693};
  auto task923 = make_shared<Task923>(tensor923, pindex);
  task916->add_dep(task923);
  task923->add_dep(task108);
  residualq->add_task(task923);

  auto tensor924 = vector<shared_ptr<Tensor>>{I1693, t2};
  auto task924 = make_shared<Task924>(tensor924, pindex);
  task923->add_dep(task924);
  task924->add_dep(task108);
  residualq->add_task(task924);

  vector<IndexRange> I1697_index = {closed_, virt_, closed_, virt_};
  auto I1697 = make_shared<Tensor>(I1697_index);
  auto tensor925 = vector<shared_ptr<Tensor>>{I1112, Gamma560_(), I1697};
  auto task925 = make_shared<Task925>(tensor925, pindex);
  task916->add_dep(task925);
  task925->add_dep(task108);
  residualq->add_task(task925);

  auto tensor926 = vector<shared_ptr<Tensor>>{I1697, t2};
  auto task926 = make_shared<Task926>(tensor926, pindex);
  task925->add_dep(task926);
  task926->add_dep(task108);
  residualq->add_task(task926);

  vector<IndexRange> I1636_index = {active_, active_, virt_, virt_};
  auto I1636 = make_shared<Tensor>(I1636_index);
  auto tensor927 = vector<shared_ptr<Tensor>>{r, I1636};
  auto task927 = make_shared<Task927>(tensor927, pindex);
  task927->add_dep(task108);
  residualq->add_task(task927);

  vector<IndexRange> I1637_index = {closed_, closed_, active_, active_};
  auto I1637 = make_shared<Tensor>(I1637_index);
  auto tensor928 = vector<shared_ptr<Tensor>>{I1636, t2, I1637};
  auto task928 = make_shared<Task928>(tensor928, pindex);
  task927->add_dep(task928);
  task928->add_dep(task108);
  residualq->add_task(task928);

  auto tensor929 = vector<shared_ptr<Tensor>>{I1637, Gamma503_(), v2_};
  auto task929 = make_shared<Task929>(tensor929, pindex);
  task928->add_dep(task929);
  task929->add_dep(task108);
  residualq->add_task(task929);

  vector<IndexRange> I1670_index = {virt_, virt_, active_, active_};
  auto I1670 = make_shared<Tensor>(I1670_index);
  auto tensor930 = vector<shared_ptr<Tensor>>{I1636, Gamma503_(), I1670};
  auto task930 = make_shared<Task930>(tensor930, pindex);
  task927->add_dep(task930);
  task930->add_dep(task108);
  residualq->add_task(task930);

  auto tensor931 = vector<shared_ptr<Tensor>>{I1670, t2, v2_};
  auto task931 = make_shared<Task931>(tensor931, pindex);
  task930->add_dep(task931);
  task931->add_dep(task108);
  residualq->add_task(task931);

  auto tensor932 = vector<shared_ptr<Tensor>>{I1636, Gamma566_(), t2};
  auto task932 = make_shared<Task932>(tensor932, pindex);
  task927->add_dep(task932);
  task932->add_dep(task108);
  residualq->add_task(task932);

  auto tensor933 = vector<shared_ptr<Tensor>>{I1636, Gamma567_(), t2};
  auto task933 = make_shared<Task933>(tensor933, pindex);
  task927->add_dep(task933);
  task933->add_dep(task108);
  residualq->add_task(task933);
}

#endif
