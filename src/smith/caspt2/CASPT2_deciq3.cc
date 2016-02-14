//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_deciq3.cc
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks15.h>
#include <src/smith/caspt2/CASPT2_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::make_deciq3(shared_ptr<Queue> deciq, shared_ptr<Task> task538, shared_ptr<Task> task539, const bool diagonal, shared_ptr<Tensor> I698) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I895_index = {active_, active_, active_, active_, active_, active_};
  auto I895 = make_shared<Tensor>(I895_index);
  auto tensor738 = vector<shared_ptr<Tensor>>{I698, Gamma299_(), I895};
  auto task738 = make_shared<Task738>(tensor738, cindex);
  task539->add_dep(task738);
  task738->add_dep(task538);
  deciq->add_task(task738);

  vector<IndexRange> I896_index = {active_, active_, virt_, active_};
  auto I896 = make_shared<Tensor>(I896_index);
  auto tensor739 = vector<shared_ptr<Tensor>>{I895, t2, I896};
  auto task739 = make_shared<Task739>(tensor739, cindex);
  task738->add_dep(task739);
  task739->add_dep(task538);
  deciq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I896, f1_, t2};
  auto task740 = make_shared<Task740>(tensor740, cindex);
  task739->add_dep(task740);
  task740->add_dep(task538);
  deciq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I895, v2_, t2};
  auto task741 = make_shared<Task741>(tensor741, cindex);
  task738->add_dep(task741);
  task741->add_dep(task538);
  deciq->add_task(task741);

  vector<IndexRange> I915_index = {active_, active_, active_, active_, active_, active_};
  auto I915 = make_shared<Tensor>(I915_index);
  auto tensor742 = vector<shared_ptr<Tensor>>{I698, Gamma304_(), I915};
  auto task742 = make_shared<Task742>(tensor742, cindex);
  task539->add_dep(task742);
  task742->add_dep(task538);
  deciq->add_task(task742);

  vector<IndexRange> I916_index = {active_, active_, virt_, active_};
  auto I916 = make_shared<Tensor>(I916_index);
  auto tensor743 = vector<shared_ptr<Tensor>>{I915, t2, I916};
  auto task743 = make_shared<Task743>(tensor743, cindex);
  task742->add_dep(task743);
  task743->add_dep(task538);
  deciq->add_task(task743);

  auto tensor744 = vector<shared_ptr<Tensor>>{I916, t2, f1_};
  auto task744 = make_shared<Task744>(tensor744, cindex);
  task743->add_dep(task744);
  task744->add_dep(task538);
  deciq->add_task(task744);

  vector<IndexRange> I919_index = {active_, active_, active_, active_, active_, active_};
  auto I919 = make_shared<Tensor>(I919_index);
  auto tensor745 = vector<shared_ptr<Tensor>>{I698, Gamma305_(), I919};
  auto task745 = make_shared<Task745>(tensor745, cindex);
  task539->add_dep(task745);
  task745->add_dep(task538);
  deciq->add_task(task745);

  vector<IndexRange> I920_index = {active_, virt_, active_, active_};
  auto I920 = make_shared<Tensor>(I920_index);
  auto tensor746 = vector<shared_ptr<Tensor>>{I919, t2, I920};
  auto task746 = make_shared<Task746>(tensor746, cindex);
  task745->add_dep(task746);
  task746->add_dep(task538);
  deciq->add_task(task746);

  auto tensor747 = vector<shared_ptr<Tensor>>{I920, t2, f1_};
  auto task747 = make_shared<Task747>(tensor747, cindex);
  task746->add_dep(task747);
  task747->add_dep(task538);
  deciq->add_task(task747);

  auto tensor748 = vector<shared_ptr<Tensor>>{I919, v2_, t2};
  auto task748 = make_shared<Task748>(tensor748, cindex);
  task745->add_dep(task748);
  task748->add_dep(task538);
  deciq->add_task(task748);

  vector<IndexRange> I923_index = {active_, active_, active_, active_, active_, active_};
  auto I923 = make_shared<Tensor>(I923_index);
  auto tensor749 = vector<shared_ptr<Tensor>>{I698, Gamma306_(), I923};
  auto task749 = make_shared<Task749>(tensor749, cindex);
  task539->add_dep(task749);
  task749->add_dep(task538);
  deciq->add_task(task749);

  auto tensor750 = vector<shared_ptr<Tensor>>{I923, t2};
  auto task750 = make_shared<Task750>(tensor750, cindex);
  task749->add_dep(task750);
  task750->add_dep(task538);
  deciq->add_task(task750);

  vector<IndexRange> I926_index = {active_, active_, active_, active_, active_, active_};
  auto I926 = make_shared<Tensor>(I926_index);
  auto tensor751 = vector<shared_ptr<Tensor>>{I698, Gamma307_(), I926};
  auto task751 = make_shared<Task751>(tensor751, cindex);
  task539->add_dep(task751);
  task751->add_dep(task538);
  deciq->add_task(task751);

  vector<IndexRange> I927_index = {active_, active_, active_, virt_};
  auto I927 = make_shared<Tensor>(I927_index);
  auto tensor752 = vector<shared_ptr<Tensor>>{I926, t2, I927};
  auto task752 = make_shared<Task752>(tensor752, cindex);
  task751->add_dep(task752);
  task752->add_dep(task538);
  deciq->add_task(task752);

  auto tensor753 = vector<shared_ptr<Tensor>>{I927, f1_, t2};
  auto task753 = make_shared<Task753>(tensor753, cindex);
  task752->add_dep(task753);
  task753->add_dep(task538);
  deciq->add_task(task753);

  vector<IndexRange> I939_index = {active_, active_, virt_, active_};
  auto I939 = make_shared<Tensor>(I939_index);
  auto tensor754 = vector<shared_ptr<Tensor>>{I926, t2, I939};
  auto task754 = make_shared<Task754>(tensor754, cindex);
  task751->add_dep(task754);
  task754->add_dep(task538);
  deciq->add_task(task754);

  auto tensor755 = vector<shared_ptr<Tensor>>{I939, t2, f1_};
  auto task755 = make_shared<Task755>(tensor755, cindex);
  task754->add_dep(task755);
  task755->add_dep(task538);
  deciq->add_task(task755);

  vector<IndexRange> I1047_index = {active_, virt_, active_, active_};
  auto I1047 = make_shared<Tensor>(I1047_index);
  auto tensor756 = vector<shared_ptr<Tensor>>{I926, t2, I1047};
  auto task756 = make_shared<Task756>(tensor756, cindex);
  task751->add_dep(task756);
  task756->add_dep(task538);
  deciq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I1047, t2};
  auto task757 = make_shared<Task757>(tensor757, cindex, this->e0_);
  task756->add_dep(task757);
  task757->add_dep(task538);
  deciq->add_task(task757);

  auto tensor758 = vector<shared_ptr<Tensor>>{I1047, f1_, t2};
  auto task758 = make_shared<Task758>(tensor758, cindex);
  task756->add_dep(task758);
  task758->add_dep(task538);
  deciq->add_task(task758);

  auto tensor759 = vector<shared_ptr<Tensor>>{I926, v2_, t2};
  auto task759 = make_shared<Task759>(tensor759, cindex);
  task751->add_dep(task759);
  task759->add_dep(task538);
  deciq->add_task(task759);

  auto tensor760 = vector<shared_ptr<Tensor>>{I926, v2_, t2};
  auto task760 = make_shared<Task760>(tensor760, cindex);
  task751->add_dep(task760);
  task760->add_dep(task538);
  deciq->add_task(task760);

  vector<IndexRange> I930_index = {active_, active_, active_, active_};
  auto I930 = make_shared<Tensor>(I930_index);
  auto tensor761 = vector<shared_ptr<Tensor>>{I698, Gamma308_(), I930};
  auto task761 = make_shared<Task761>(tensor761, cindex);
  task539->add_dep(task761);
  task761->add_dep(task538);
  deciq->add_task(task761);

  vector<IndexRange> I931_index = {active_, virt_};
  auto I931 = make_shared<Tensor>(I931_index);
  auto tensor762 = vector<shared_ptr<Tensor>>{I930, t2, I931};
  auto task762 = make_shared<Task762>(tensor762, cindex);
  task761->add_dep(task762);
  task762->add_dep(task538);
  deciq->add_task(task762);

  auto tensor763 = vector<shared_ptr<Tensor>>{I931, t2, f1_};
  auto task763 = make_shared<Task763>(tensor763, cindex);
  task762->add_dep(task763);
  task763->add_dep(task538);
  deciq->add_task(task763);

  auto tensor764 = vector<shared_ptr<Tensor>>{I931, t2, f1_};
  auto task764 = make_shared<Task764>(tensor764, cindex);
  task762->add_dep(task764);
  task764->add_dep(task538);
  deciq->add_task(task764);

  vector<IndexRange> I997_index = {virt_, active_};
  auto I997 = make_shared<Tensor>(I997_index);
  auto tensor765 = vector<shared_ptr<Tensor>>{I930, t2, I997};
  auto task765 = make_shared<Task765>(tensor765, cindex);
  task761->add_dep(task765);
  task765->add_dep(task538);
  deciq->add_task(task765);

  auto tensor766 = vector<shared_ptr<Tensor>>{I997, f1_, t2};
  auto task766 = make_shared<Task766>(tensor766, cindex);
  task765->add_dep(task766);
  task766->add_dep(task538);
  deciq->add_task(task766);

  vector<IndexRange> I1001_index = {virt_, active_};
  auto I1001 = make_shared<Tensor>(I1001_index);
  auto tensor767 = vector<shared_ptr<Tensor>>{I930, t2, I1001};
  auto task767 = make_shared<Task767>(tensor767, cindex);
  task761->add_dep(task767);
  task767->add_dep(task538);
  deciq->add_task(task767);

  auto tensor768 = vector<shared_ptr<Tensor>>{I1001, f1_, t2};
  auto task768 = make_shared<Task768>(tensor768, cindex);
  task767->add_dep(task768);
  task768->add_dep(task538);
  deciq->add_task(task768);

  vector<IndexRange> I1043_index = {virt_, virt_, active_, active_};
  auto I1043 = make_shared<Tensor>(I1043_index);
  auto tensor769 = vector<shared_ptr<Tensor>>{I930, t2, I1043};
  auto task769 = make_shared<Task769>(tensor769, cindex);
  task761->add_dep(task769);
  task769->add_dep(task538);
  deciq->add_task(task769);

  auto tensor770 = vector<shared_ptr<Tensor>>{I1043, f1_, t2};
  auto task770 = make_shared<Task770>(tensor770, cindex);
  task769->add_dep(task770);
  task770->add_dep(task538);
  deciq->add_task(task770);

  auto tensor771 = vector<shared_ptr<Tensor>>{I1043, f1_, t2};
  auto task771 = make_shared<Task771>(tensor771, cindex);
  task769->add_dep(task771);
  task771->add_dep(task538);
  deciq->add_task(task771);

  vector<IndexRange> I1051_index = {active_, active_, virt_, virt_};
  auto I1051 = make_shared<Tensor>(I1051_index);
  auto tensor772 = vector<shared_ptr<Tensor>>{I930, t2, I1051};
  auto task772 = make_shared<Task772>(tensor772, cindex);
  task761->add_dep(task772);
  task772->add_dep(task538);
  deciq->add_task(task772);

  auto tensor773 = vector<shared_ptr<Tensor>>{I1051, t2, f1_};
  auto task773 = make_shared<Task773>(tensor773, cindex);
  task772->add_dep(task773);
  task773->add_dep(task538);
  deciq->add_task(task773);

  auto tensor774 = vector<shared_ptr<Tensor>>{I930, t2};
  auto task774 = make_shared<Task774>(tensor774, cindex, this->e0_);
  task761->add_dep(task774);
  task774->add_dep(task538);
  deciq->add_task(task774);

  auto tensor775 = vector<shared_ptr<Tensor>>{I930, v2_, t2};
  auto task775 = make_shared<Task775>(tensor775, cindex);
  task761->add_dep(task775);
  task775->add_dep(task538);
  deciq->add_task(task775);

  auto tensor776 = vector<shared_ptr<Tensor>>{I930, v2_, t2};
  auto task776 = make_shared<Task776>(tensor776, cindex);
  task761->add_dep(task776);
  task776->add_dep(task538);
  deciq->add_task(task776);

  auto tensor777 = vector<shared_ptr<Tensor>>{I930, h1_, t2};
  auto task777 = make_shared<Task777>(tensor777, cindex);
  task761->add_dep(task777);
  task777->add_dep(task538);
  deciq->add_task(task777);

  auto tensor778 = vector<shared_ptr<Tensor>>{I930, h1_, t2};
  auto task778 = make_shared<Task778>(tensor778, cindex);
  task761->add_dep(task778);
  task778->add_dep(task538);
  deciq->add_task(task778);

  vector<IndexRange> I966_index;
  auto I966 = make_shared<Tensor>(I966_index);
  auto tensor779 = vector<shared_ptr<Tensor>>{I698, Gamma317_(), I966};
  auto task779 = make_shared<Task779>(tensor779, cindex);
  task539->add_dep(task779);
  task779->add_dep(task538);
  deciq->add_task(task779);

  auto tensor780 = vector<shared_ptr<Tensor>>{I966, t2};
  auto task780 = make_shared<Task780>(tensor780, cindex);
  task779->add_dep(task780);
  task780->add_dep(task538);
  deciq->add_task(task780);

  auto tensor781 = vector<shared_ptr<Tensor>>{I966, t2};
  auto task781 = make_shared<Task781>(tensor781, cindex);
  task779->add_dep(task781);
  task781->add_dep(task538);
  deciq->add_task(task781);

  vector<IndexRange> I1012_index = {active_, active_};
  auto I1012 = make_shared<Tensor>(I1012_index);
  auto tensor782 = vector<shared_ptr<Tensor>>{I698, Gamma329_(), I1012};
  auto task782 = make_shared<Task782>(tensor782, cindex);
  task539->add_dep(task782);
  task782->add_dep(task538);
  deciq->add_task(task782);

  auto tensor783 = vector<shared_ptr<Tensor>>{I1012, t2};
  auto task783 = make_shared<Task783>(tensor783, cindex);
  task782->add_dep(task783);
  task783->add_dep(task538);
  deciq->add_task(task783);

  auto tensor784 = vector<shared_ptr<Tensor>>{I1012, t2};
  auto task784 = make_shared<Task784>(tensor784, cindex);
  task782->add_dep(task784);
  task784->add_dep(task538);
  deciq->add_task(task784);

  vector<IndexRange> I1054_index = {active_, active_, active_, active_};
  auto I1054 = make_shared<Tensor>(I1054_index);
  auto tensor785 = vector<shared_ptr<Tensor>>{I698, Gamma340_(), I1054};
  auto task785 = make_shared<Task785>(tensor785, cindex);
  task539->add_dep(task785);
  task785->add_dep(task538);
  deciq->add_task(task785);

  auto tensor786 = vector<shared_ptr<Tensor>>{I1054, t2};
  auto task786 = make_shared<Task786>(tensor786, cindex);
  task785->add_dep(task786);
  task786->add_dep(task538);
  deciq->add_task(task786);

  vector<IndexRange> I1100_index = {active_, active_, active_, active_, active_, active_};
  auto I1100 = make_shared<Tensor>(I1100_index);
  auto tensor787 = vector<shared_ptr<Tensor>>{I698, Gamma355_(), I1100};
  auto task787 = make_shared<Task787>(tensor787, cindex);
  task539->add_dep(task787);
  task787->add_dep(task538);
  deciq->add_task(task787);

  auto tensor788 = vector<shared_ptr<Tensor>>{I1100, v2_, t2};
  auto task788 = make_shared<Task788>(tensor788, cindex);
  task787->add_dep(task788);
  task788->add_dep(task538);
  deciq->add_task(task788);

  vector<IndexRange> I1154_index = {active_, active_, active_, active_, active_, active_};
  auto I1154 = make_shared<Tensor>(I1154_index);
  auto tensor789 = vector<shared_ptr<Tensor>>{I698, Gamma373_(), I1154};
  auto task789 = make_shared<Task789>(tensor789, cindex);
  task539->add_dep(task789);
  task789->add_dep(task538);
  deciq->add_task(task789);

  auto tensor790 = vector<shared_ptr<Tensor>>{I1154, v2_, t2};
  auto task790 = make_shared<Task790>(tensor790, cindex);
  task789->add_dep(task790);
  task790->add_dep(task538);
  deciq->add_task(task790);
}

#endif
