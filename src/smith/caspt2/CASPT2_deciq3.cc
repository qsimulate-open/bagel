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

void CASPT2::CASPT2::make_deciq3(shared_ptr<Queue> deciq, shared_ptr<Task> task539, shared_ptr<Task> task540, const bool diagonal, shared_ptr<Tensor> I684) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I881_index = {active_, active_, active_, active_, active_, active_};
  auto I881 = make_shared<Tensor>(I881_index);
  auto tensor739 = vector<shared_ptr<Tensor>>{I684, Gamma299_(), I881};
  auto task739 = make_shared<Task739>(tensor739, cindex);
  task540->add_dep(task739);
  task739->add_dep(task539);
  deciq->add_task(task739);

  vector<IndexRange> I882_index = {active_, active_, virt_, active_};
  auto I882 = make_shared<Tensor>(I882_index);
  auto tensor740 = vector<shared_ptr<Tensor>>{I881, t2, I882};
  auto task740 = make_shared<Task740>(tensor740, cindex);
  task739->add_dep(task740);
  task740->add_dep(task539);
  deciq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I882, f1_, t2};
  auto task741 = make_shared<Task741>(tensor741, cindex);
  task740->add_dep(task741);
  task741->add_dep(task539);
  deciq->add_task(task741);

  auto tensor742 = vector<shared_ptr<Tensor>>{I881, v2_, t2};
  auto task742 = make_shared<Task742>(tensor742, cindex);
  task739->add_dep(task742);
  task742->add_dep(task539);
  deciq->add_task(task742);

  vector<IndexRange> I901_index = {active_, active_, active_, active_, active_, active_};
  auto I901 = make_shared<Tensor>(I901_index);
  auto tensor743 = vector<shared_ptr<Tensor>>{I684, Gamma304_(), I901};
  auto task743 = make_shared<Task743>(tensor743, cindex);
  task540->add_dep(task743);
  task743->add_dep(task539);
  deciq->add_task(task743);

  vector<IndexRange> I902_index = {active_, active_, virt_, active_};
  auto I902 = make_shared<Tensor>(I902_index);
  auto tensor744 = vector<shared_ptr<Tensor>>{I901, t2, I902};
  auto task744 = make_shared<Task744>(tensor744, cindex);
  task743->add_dep(task744);
  task744->add_dep(task539);
  deciq->add_task(task744);

  auto tensor745 = vector<shared_ptr<Tensor>>{I902, t2, f1_};
  auto task745 = make_shared<Task745>(tensor745, cindex);
  task744->add_dep(task745);
  task745->add_dep(task539);
  deciq->add_task(task745);

  vector<IndexRange> I905_index = {active_, active_, active_, active_, active_, active_};
  auto I905 = make_shared<Tensor>(I905_index);
  auto tensor746 = vector<shared_ptr<Tensor>>{I684, Gamma305_(), I905};
  auto task746 = make_shared<Task746>(tensor746, cindex);
  task540->add_dep(task746);
  task746->add_dep(task539);
  deciq->add_task(task746);

  vector<IndexRange> I906_index = {active_, virt_, active_, active_};
  auto I906 = make_shared<Tensor>(I906_index);
  auto tensor747 = vector<shared_ptr<Tensor>>{I905, t2, I906};
  auto task747 = make_shared<Task747>(tensor747, cindex);
  task746->add_dep(task747);
  task747->add_dep(task539);
  deciq->add_task(task747);

  auto tensor748 = vector<shared_ptr<Tensor>>{I906, t2, f1_};
  auto task748 = make_shared<Task748>(tensor748, cindex);
  task747->add_dep(task748);
  task748->add_dep(task539);
  deciq->add_task(task748);

  auto tensor749 = vector<shared_ptr<Tensor>>{I905, v2_, t2};
  auto task749 = make_shared<Task749>(tensor749, cindex);
  task746->add_dep(task749);
  task749->add_dep(task539);
  deciq->add_task(task749);

  vector<IndexRange> I909_index = {active_, active_, active_, active_, active_, active_};
  auto I909 = make_shared<Tensor>(I909_index);
  auto tensor750 = vector<shared_ptr<Tensor>>{I684, Gamma306_(), I909};
  auto task750 = make_shared<Task750>(tensor750, cindex);
  task540->add_dep(task750);
  task750->add_dep(task539);
  deciq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I909, t2};
  auto task751 = make_shared<Task751>(tensor751, cindex);
  task750->add_dep(task751);
  task751->add_dep(task539);
  deciq->add_task(task751);

  vector<IndexRange> I912_index = {active_, active_, active_, active_, active_, active_};
  auto I912 = make_shared<Tensor>(I912_index);
  auto tensor752 = vector<shared_ptr<Tensor>>{I684, Gamma307_(), I912};
  auto task752 = make_shared<Task752>(tensor752, cindex);
  task540->add_dep(task752);
  task752->add_dep(task539);
  deciq->add_task(task752);

  vector<IndexRange> I913_index = {active_, active_, active_, virt_};
  auto I913 = make_shared<Tensor>(I913_index);
  auto tensor753 = vector<shared_ptr<Tensor>>{I912, t2, I913};
  auto task753 = make_shared<Task753>(tensor753, cindex);
  task752->add_dep(task753);
  task753->add_dep(task539);
  deciq->add_task(task753);

  auto tensor754 = vector<shared_ptr<Tensor>>{I913, f1_, t2};
  auto task754 = make_shared<Task754>(tensor754, cindex);
  task753->add_dep(task754);
  task754->add_dep(task539);
  deciq->add_task(task754);

  vector<IndexRange> I925_index = {active_, active_, virt_, active_};
  auto I925 = make_shared<Tensor>(I925_index);
  auto tensor755 = vector<shared_ptr<Tensor>>{I912, t2, I925};
  auto task755 = make_shared<Task755>(tensor755, cindex);
  task752->add_dep(task755);
  task755->add_dep(task539);
  deciq->add_task(task755);

  auto tensor756 = vector<shared_ptr<Tensor>>{I925, t2, f1_};
  auto task756 = make_shared<Task756>(tensor756, cindex);
  task755->add_dep(task756);
  task756->add_dep(task539);
  deciq->add_task(task756);

  vector<IndexRange> I1033_index = {active_, virt_, active_, active_};
  auto I1033 = make_shared<Tensor>(I1033_index);
  auto tensor757 = vector<shared_ptr<Tensor>>{I912, t2, I1033};
  auto task757 = make_shared<Task757>(tensor757, cindex);
  task752->add_dep(task757);
  task757->add_dep(task539);
  deciq->add_task(task757);

  auto tensor758 = vector<shared_ptr<Tensor>>{I1033, t2};
  auto task758 = make_shared<Task758>(tensor758, cindex, this->e0_);
  task757->add_dep(task758);
  task758->add_dep(task539);
  deciq->add_task(task758);

  auto tensor759 = vector<shared_ptr<Tensor>>{I1033, f1_, t2};
  auto task759 = make_shared<Task759>(tensor759, cindex);
  task757->add_dep(task759);
  task759->add_dep(task539);
  deciq->add_task(task759);

  auto tensor760 = vector<shared_ptr<Tensor>>{I912, v2_, t2};
  auto task760 = make_shared<Task760>(tensor760, cindex);
  task752->add_dep(task760);
  task760->add_dep(task539);
  deciq->add_task(task760);

  auto tensor761 = vector<shared_ptr<Tensor>>{I912, v2_, t2};
  auto task761 = make_shared<Task761>(tensor761, cindex);
  task752->add_dep(task761);
  task761->add_dep(task539);
  deciq->add_task(task761);

  vector<IndexRange> I916_index = {active_, active_, active_, active_};
  auto I916 = make_shared<Tensor>(I916_index);
  auto tensor762 = vector<shared_ptr<Tensor>>{I684, Gamma308_(), I916};
  auto task762 = make_shared<Task762>(tensor762, cindex);
  task540->add_dep(task762);
  task762->add_dep(task539);
  deciq->add_task(task762);

  vector<IndexRange> I917_index = {active_, virt_};
  auto I917 = make_shared<Tensor>(I917_index);
  auto tensor763 = vector<shared_ptr<Tensor>>{I916, t2, I917};
  auto task763 = make_shared<Task763>(tensor763, cindex);
  task762->add_dep(task763);
  task763->add_dep(task539);
  deciq->add_task(task763);

  auto tensor764 = vector<shared_ptr<Tensor>>{I917, t2, f1_};
  auto task764 = make_shared<Task764>(tensor764, cindex);
  task763->add_dep(task764);
  task764->add_dep(task539);
  deciq->add_task(task764);

  auto tensor765 = vector<shared_ptr<Tensor>>{I917, t2, f1_};
  auto task765 = make_shared<Task765>(tensor765, cindex);
  task763->add_dep(task765);
  task765->add_dep(task539);
  deciq->add_task(task765);

  vector<IndexRange> I983_index = {virt_, active_};
  auto I983 = make_shared<Tensor>(I983_index);
  auto tensor766 = vector<shared_ptr<Tensor>>{I916, t2, I983};
  auto task766 = make_shared<Task766>(tensor766, cindex);
  task762->add_dep(task766);
  task766->add_dep(task539);
  deciq->add_task(task766);

  auto tensor767 = vector<shared_ptr<Tensor>>{I983, f1_, t2};
  auto task767 = make_shared<Task767>(tensor767, cindex);
  task766->add_dep(task767);
  task767->add_dep(task539);
  deciq->add_task(task767);

  vector<IndexRange> I987_index = {virt_, active_};
  auto I987 = make_shared<Tensor>(I987_index);
  auto tensor768 = vector<shared_ptr<Tensor>>{I916, t2, I987};
  auto task768 = make_shared<Task768>(tensor768, cindex);
  task762->add_dep(task768);
  task768->add_dep(task539);
  deciq->add_task(task768);

  auto tensor769 = vector<shared_ptr<Tensor>>{I987, f1_, t2};
  auto task769 = make_shared<Task769>(tensor769, cindex);
  task768->add_dep(task769);
  task769->add_dep(task539);
  deciq->add_task(task769);

  vector<IndexRange> I1029_index = {virt_, virt_, active_, active_};
  auto I1029 = make_shared<Tensor>(I1029_index);
  auto tensor770 = vector<shared_ptr<Tensor>>{I916, t2, I1029};
  auto task770 = make_shared<Task770>(tensor770, cindex);
  task762->add_dep(task770);
  task770->add_dep(task539);
  deciq->add_task(task770);

  auto tensor771 = vector<shared_ptr<Tensor>>{I1029, f1_, t2};
  auto task771 = make_shared<Task771>(tensor771, cindex);
  task770->add_dep(task771);
  task771->add_dep(task539);
  deciq->add_task(task771);

  auto tensor772 = vector<shared_ptr<Tensor>>{I1029, f1_, t2};
  auto task772 = make_shared<Task772>(tensor772, cindex);
  task770->add_dep(task772);
  task772->add_dep(task539);
  deciq->add_task(task772);

  vector<IndexRange> I1037_index = {active_, active_, virt_, virt_};
  auto I1037 = make_shared<Tensor>(I1037_index);
  auto tensor773 = vector<shared_ptr<Tensor>>{I916, t2, I1037};
  auto task773 = make_shared<Task773>(tensor773, cindex);
  task762->add_dep(task773);
  task773->add_dep(task539);
  deciq->add_task(task773);

  auto tensor774 = vector<shared_ptr<Tensor>>{I1037, t2, f1_};
  auto task774 = make_shared<Task774>(tensor774, cindex);
  task773->add_dep(task774);
  task774->add_dep(task539);
  deciq->add_task(task774);

  auto tensor775 = vector<shared_ptr<Tensor>>{I916, t2};
  auto task775 = make_shared<Task775>(tensor775, cindex, this->e0_);
  task762->add_dep(task775);
  task775->add_dep(task539);
  deciq->add_task(task775);

  auto tensor776 = vector<shared_ptr<Tensor>>{I916, v2_, t2};
  auto task776 = make_shared<Task776>(tensor776, cindex);
  task762->add_dep(task776);
  task776->add_dep(task539);
  deciq->add_task(task776);

  auto tensor777 = vector<shared_ptr<Tensor>>{I916, v2_, t2};
  auto task777 = make_shared<Task777>(tensor777, cindex);
  task762->add_dep(task777);
  task777->add_dep(task539);
  deciq->add_task(task777);

  auto tensor778 = vector<shared_ptr<Tensor>>{I916, h1_, t2};
  auto task778 = make_shared<Task778>(tensor778, cindex);
  task762->add_dep(task778);
  task778->add_dep(task539);
  deciq->add_task(task778);

  auto tensor779 = vector<shared_ptr<Tensor>>{I916, h1_, t2};
  auto task779 = make_shared<Task779>(tensor779, cindex);
  task762->add_dep(task779);
  task779->add_dep(task539);
  deciq->add_task(task779);

  vector<IndexRange> I952_index;
  auto I952 = make_shared<Tensor>(I952_index);
  auto tensor780 = vector<shared_ptr<Tensor>>{I684, Gamma317_(), I952};
  auto task780 = make_shared<Task780>(tensor780, cindex);
  task540->add_dep(task780);
  task780->add_dep(task539);
  deciq->add_task(task780);

  auto tensor781 = vector<shared_ptr<Tensor>>{I952, t2};
  auto task781 = make_shared<Task781>(tensor781, cindex);
  task780->add_dep(task781);
  task781->add_dep(task539);
  deciq->add_task(task781);

  auto tensor782 = vector<shared_ptr<Tensor>>{I952, t2};
  auto task782 = make_shared<Task782>(tensor782, cindex);
  task780->add_dep(task782);
  task782->add_dep(task539);
  deciq->add_task(task782);

  vector<IndexRange> I998_index = {active_, active_};
  auto I998 = make_shared<Tensor>(I998_index);
  auto tensor783 = vector<shared_ptr<Tensor>>{I684, Gamma329_(), I998};
  auto task783 = make_shared<Task783>(tensor783, cindex);
  task540->add_dep(task783);
  task783->add_dep(task539);
  deciq->add_task(task783);

  auto tensor784 = vector<shared_ptr<Tensor>>{I998, t2};
  auto task784 = make_shared<Task784>(tensor784, cindex);
  task783->add_dep(task784);
  task784->add_dep(task539);
  deciq->add_task(task784);

  auto tensor785 = vector<shared_ptr<Tensor>>{I998, t2};
  auto task785 = make_shared<Task785>(tensor785, cindex);
  task783->add_dep(task785);
  task785->add_dep(task539);
  deciq->add_task(task785);

  vector<IndexRange> I1040_index = {active_, active_, active_, active_};
  auto I1040 = make_shared<Tensor>(I1040_index);
  auto tensor786 = vector<shared_ptr<Tensor>>{I684, Gamma340_(), I1040};
  auto task786 = make_shared<Task786>(tensor786, cindex);
  task540->add_dep(task786);
  task786->add_dep(task539);
  deciq->add_task(task786);

  auto tensor787 = vector<shared_ptr<Tensor>>{I1040, t2};
  auto task787 = make_shared<Task787>(tensor787, cindex);
  task786->add_dep(task787);
  task787->add_dep(task539);
  deciq->add_task(task787);

  vector<IndexRange> I1086_index = {active_, active_, active_, active_, active_, active_};
  auto I1086 = make_shared<Tensor>(I1086_index);
  auto tensor788 = vector<shared_ptr<Tensor>>{I684, Gamma355_(), I1086};
  auto task788 = make_shared<Task788>(tensor788, cindex);
  task540->add_dep(task788);
  task788->add_dep(task539);
  deciq->add_task(task788);

  auto tensor789 = vector<shared_ptr<Tensor>>{I1086, v2_, t2};
  auto task789 = make_shared<Task789>(tensor789, cindex);
  task788->add_dep(task789);
  task789->add_dep(task539);
  deciq->add_task(task789);

  vector<IndexRange> I1140_index = {active_, active_, active_, active_, active_, active_};
  auto I1140 = make_shared<Tensor>(I1140_index);
  auto tensor790 = vector<shared_ptr<Tensor>>{I684, Gamma373_(), I1140};
  auto task790 = make_shared<Task790>(tensor790, cindex);
  task540->add_dep(task790);
  task790->add_dep(task539);
  deciq->add_task(task790);

  auto tensor791 = vector<shared_ptr<Tensor>>{I1140, v2_, t2};
  auto task791 = make_shared<Task791>(tensor791, cindex);
  task790->add_dep(task791);
  task791->add_dep(task539);
  deciq->add_task(task791);
}

#endif
