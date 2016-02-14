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
#include <src/smith/mrci/MRCI_tasks15.h>
#include <src/smith/mrci/MRCI_tasks16.h>
#include <src/smith/mrci/MRCI_tasks17.h>
#include <src/smith/mrci/MRCI_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq8(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I194_index = {virt_, closed_, active_, virt_};
  auto I194 = make_shared<Tensor>(I194_index);
  auto tensor727 = vector<shared_ptr<Tensor>>{r, I194};
  auto task727 = make_shared<Task727>(tensor727, pindex);
  task727->add_dep(task108);
  residualq->add_task(task727);

  vector<IndexRange> I195_index = {virt_, closed_, active_, active_};
  auto I195 = make_shared<Tensor>(I195_index);
  auto tensor728 = vector<shared_ptr<Tensor>>{I194, h1_, I195};
  auto task728 = make_shared<Task728>(tensor728, pindex);
  task727->add_dep(task728);
  task728->add_dep(task108);
  residualq->add_task(task728);

  vector<IndexRange> I196_index = {active_, virt_, closed_, active_};
  auto I196 = make_shared<Tensor>(I196_index);
  auto tensor729 = vector<shared_ptr<Tensor>>{I195, Gamma29_(), I196};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task728->add_dep(task729);
  task729->add_dep(task108);
  residualq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I196, t2};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task729->add_dep(task730);
  task730->add_dep(task108);
  residualq->add_task(task730);

  vector<IndexRange> I198_index = {virt_, closed_, active_, active_};
  auto I198 = make_shared<Tensor>(I198_index);
  auto tensor731 = vector<shared_ptr<Tensor>>{I194, h1_, I198};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task727->add_dep(task731);
  task731->add_dep(task108);
  residualq->add_task(task731);

  auto tensor732 = vector<shared_ptr<Tensor>>{I198, Gamma27_(), t2};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task731->add_dep(task732);
  task732->add_dep(task108);
  residualq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I198, Gamma29_(), t2};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task731->add_dep(task733);
  task733->add_dep(task108);
  residualq->add_task(task733);

  vector<IndexRange> I207_index = {virt_, active_};
  auto I207 = make_shared<Tensor>(I207_index);
  auto tensor734 = vector<shared_ptr<Tensor>>{I194, h1_, I207};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task727->add_dep(task734);
  task734->add_dep(task108);
  residualq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I207, Gamma51_(), t2};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task734->add_dep(task735);
  task735->add_dep(task108);
  residualq->add_task(task735);

  vector<IndexRange> I210_index = {virt_, active_};
  auto I210 = make_shared<Tensor>(I210_index);
  auto tensor736 = vector<shared_ptr<Tensor>>{I194, h1_, I210};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task727->add_dep(task736);
  task736->add_dep(task108);
  residualq->add_task(task736);

  auto tensor737 = vector<shared_ptr<Tensor>>{I210, Gamma51_(), t2};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task736->add_dep(task737);
  task737->add_dep(task108);
  residualq->add_task(task737);

  vector<IndexRange> I213_index = {closed_, active_};
  auto I213 = make_shared<Tensor>(I213_index);
  auto tensor738 = vector<shared_ptr<Tensor>>{I194, t2, I213};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task727->add_dep(task738);
  task738->add_dep(task108);
  residualq->add_task(task738);

  auto tensor739 = vector<shared_ptr<Tensor>>{I213, Gamma32_(), h1_};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task738->add_dep(task739);
  task739->add_dep(task108);
  residualq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I213, Gamma51_(), v2_};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task738->add_dep(task740);
  task740->add_dep(task108);
  residualq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I213, Gamma29_(), v2_};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task738->add_dep(task741);
  task741->add_dep(task108);
  residualq->add_task(task741);

  vector<IndexRange> I216_index = {closed_, active_};
  auto I216 = make_shared<Tensor>(I216_index);
  auto tensor742 = vector<shared_ptr<Tensor>>{I194, t2, I216};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task727->add_dep(task742);
  task742->add_dep(task108);
  residualq->add_task(task742);

  auto tensor743 = vector<shared_ptr<Tensor>>{I216, Gamma32_(), h1_};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task108);
  residualq->add_task(task743);

  auto tensor744 = vector<shared_ptr<Tensor>>{I216, Gamma51_(), v2_};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task742->add_dep(task744);
  task744->add_dep(task108);
  residualq->add_task(task744);

  auto tensor745 = vector<shared_ptr<Tensor>>{I216, Gamma29_(), v2_};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task742->add_dep(task745);
  task745->add_dep(task108);
  residualq->add_task(task745);

  vector<IndexRange> I219_index = {closed_, active_, virt_, virt_};
  auto I219 = make_shared<Tensor>(I219_index);
  auto tensor746 = vector<shared_ptr<Tensor>>{I194, Gamma32_(), I219};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task727->add_dep(task746);
  task746->add_dep(task108);
  residualq->add_task(task746);

  auto tensor747 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task746->add_dep(task747);
  task747->add_dep(task108);
  residualq->add_task(task747);

  auto tensor748 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task748 = make_shared<Task748>(tensor748, pindex);
  task746->add_dep(task748);
  task748->add_dep(task108);
  residualq->add_task(task748);

  auto tensor749 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task746->add_dep(task749);
  task749->add_dep(task108);
  residualq->add_task(task749);

  auto tensor750 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task746->add_dep(task750);
  task750->add_dep(task108);
  residualq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task746->add_dep(task751);
  task751->add_dep(task108);
  residualq->add_task(task751);

  auto tensor752 = vector<shared_ptr<Tensor>>{I219, t2, h1_};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task746->add_dep(task752);
  task752->add_dep(task108);
  residualq->add_task(task752);

  auto tensor753 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task753 = make_shared<Task753>(tensor753, pindex);
  task746->add_dep(task753);
  task753->add_dep(task108);
  residualq->add_task(task753);

  auto tensor754 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task746->add_dep(task754);
  task754->add_dep(task108);
  residualq->add_task(task754);

  auto tensor755 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task755 = make_shared<Task755>(tensor755, pindex);
  task746->add_dep(task755);
  task755->add_dep(task108);
  residualq->add_task(task755);

  auto tensor756 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task746->add_dep(task756);
  task756->add_dep(task108);
  residualq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task746->add_dep(task757);
  task757->add_dep(task108);
  residualq->add_task(task757);

  auto tensor758 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task758 = make_shared<Task758>(tensor758, pindex);
  task746->add_dep(task758);
  task758->add_dep(task108);
  residualq->add_task(task758);

  auto tensor759 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task746->add_dep(task759);
  task759->add_dep(task108);
  residualq->add_task(task759);

  auto tensor760 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task746->add_dep(task760);
  task760->add_dep(task108);
  residualq->add_task(task760);

  auto tensor761 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task746->add_dep(task761);
  task761->add_dep(task108);
  residualq->add_task(task761);

  auto tensor762 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task762 = make_shared<Task762>(tensor762, pindex);
  task746->add_dep(task762);
  task762->add_dep(task108);
  residualq->add_task(task762);

  auto tensor763 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task763 = make_shared<Task763>(tensor763, pindex);
  task746->add_dep(task763);
  task763->add_dep(task108);
  residualq->add_task(task763);

  auto tensor764 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task764 = make_shared<Task764>(tensor764, pindex);
  task746->add_dep(task764);
  task764->add_dep(task108);
  residualq->add_task(task764);

  auto tensor765 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task765 = make_shared<Task765>(tensor765, pindex);
  task746->add_dep(task765);
  task765->add_dep(task108);
  residualq->add_task(task765);

  auto tensor766 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task766 = make_shared<Task766>(tensor766, pindex);
  task746->add_dep(task766);
  task766->add_dep(task108);
  residualq->add_task(task766);

  auto tensor767 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task767 = make_shared<Task767>(tensor767, pindex);
  task746->add_dep(task767);
  task767->add_dep(task108);
  residualq->add_task(task767);

  auto tensor768 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task768 = make_shared<Task768>(tensor768, pindex);
  task746->add_dep(task768);
  task768->add_dep(task108);
  residualq->add_task(task768);

  auto tensor769 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task769 = make_shared<Task769>(tensor769, pindex);
  task746->add_dep(task769);
  task769->add_dep(task108);
  residualq->add_task(task769);

  auto tensor770 = vector<shared_ptr<Tensor>>{I219, t2, v2_};
  auto task770 = make_shared<Task770>(tensor770, pindex);
  task746->add_dep(task770);
  task770->add_dep(task108);
  residualq->add_task(task770);

  vector<IndexRange> I237_index = {virt_, virt_, active_, active_};
  auto I237 = make_shared<Tensor>(I237_index);
  auto tensor771 = vector<shared_ptr<Tensor>>{I194, h1_, I237};
  auto task771 = make_shared<Task771>(tensor771, pindex);
  task727->add_dep(task771);
  task771->add_dep(task108);
  residualq->add_task(task771);

  auto tensor772 = vector<shared_ptr<Tensor>>{I237, Gamma51_(), t2};
  auto task772 = make_shared<Task772>(tensor772, pindex);
  task771->add_dep(task772);
  task772->add_dep(task108);
  residualq->add_task(task772);

  vector<IndexRange> I1355_index = {closed_, active_, active_, active_};
  auto I1355 = make_shared<Tensor>(I1355_index);
  auto tensor773 = vector<shared_ptr<Tensor>>{I194, v2_, I1355};
  auto task773 = make_shared<Task773>(tensor773, pindex);
  task727->add_dep(task773);
  task773->add_dep(task108);
  residualq->add_task(task773);

  auto tensor774 = vector<shared_ptr<Tensor>>{I1355, Gamma24_(), t2};
  auto task774 = make_shared<Task774>(tensor774, pindex);
  task773->add_dep(task774);
  task774->add_dep(task108);
  residualq->add_task(task774);

  vector<IndexRange> I1358_index = {virt_, closed_, active_, active_};
  auto I1358 = make_shared<Tensor>(I1358_index);
  auto tensor775 = vector<shared_ptr<Tensor>>{I194, t2, I1358};
  auto task775 = make_shared<Task775>(tensor775, pindex);
  task727->add_dep(task775);
  task775->add_dep(task108);
  residualq->add_task(task775);

  auto tensor776 = vector<shared_ptr<Tensor>>{I1358, Gamma25_(), v2_};
  auto task776 = make_shared<Task776>(tensor776, pindex);
  task775->add_dep(task776);
  task776->add_dep(task108);
  residualq->add_task(task776);

  vector<IndexRange> I1361_index = {virt_, closed_, active_, active_};
  auto I1361 = make_shared<Tensor>(I1361_index);
  auto tensor777 = vector<shared_ptr<Tensor>>{I194, t2, I1361};
  auto task777 = make_shared<Task777>(tensor777, pindex);
  task727->add_dep(task777);
  task777->add_dep(task108);
  residualq->add_task(task777);

  auto tensor778 = vector<shared_ptr<Tensor>>{I1361, Gamma5_(), v2_};
  auto task778 = make_shared<Task778>(tensor778, pindex);
  task777->add_dep(task778);
  task778->add_dep(task108);
  residualq->add_task(task778);

  vector<IndexRange> I1364_index = {virt_, closed_, active_, active_};
  auto I1364 = make_shared<Tensor>(I1364_index);
  auto tensor779 = vector<shared_ptr<Tensor>>{I194, t2, I1364};
  auto task779 = make_shared<Task779>(tensor779, pindex);
  task727->add_dep(task779);
  task779->add_dep(task108);
  residualq->add_task(task779);

  auto tensor780 = vector<shared_ptr<Tensor>>{I1364, Gamma25_(), v2_};
  auto task780 = make_shared<Task780>(tensor780, pindex);
  task779->add_dep(task780);
  task780->add_dep(task108);
  residualq->add_task(task780);

  vector<IndexRange> I1367_index = {virt_, closed_, active_, active_};
  auto I1367 = make_shared<Tensor>(I1367_index);
  auto tensor781 = vector<shared_ptr<Tensor>>{I194, t2, I1367};
  auto task781 = make_shared<Task781>(tensor781, pindex);
  task727->add_dep(task781);
  task781->add_dep(task108);
  residualq->add_task(task781);

  auto tensor782 = vector<shared_ptr<Tensor>>{I1367, Gamma25_(), v2_};
  auto task782 = make_shared<Task782>(tensor782, pindex);
  task781->add_dep(task782);
  task782->add_dep(task108);
  residualq->add_task(task782);

  vector<IndexRange> I1370_index = {virt_, active_, active_, active_};
  auto I1370 = make_shared<Tensor>(I1370_index);
  auto tensor783 = vector<shared_ptr<Tensor>>{I194, t2, I1370};
  auto task783 = make_shared<Task783>(tensor783, pindex);
  task727->add_dep(task783);
  task783->add_dep(task108);
  residualq->add_task(task783);

  auto tensor784 = vector<shared_ptr<Tensor>>{I1370, Gamma49_(), v2_};
  auto task784 = make_shared<Task784>(tensor784, pindex);
  task783->add_dep(task784);
  task784->add_dep(task108);
  residualq->add_task(task784);

  auto tensor785 = vector<shared_ptr<Tensor>>{I1370, Gamma240_(), v2_};
  auto task785 = make_shared<Task785>(tensor785, pindex);
  task783->add_dep(task785);
  task785->add_dep(task108);
  residualq->add_task(task785);

  vector<IndexRange> I1373_index = {virt_, active_, active_, active_};
  auto I1373 = make_shared<Tensor>(I1373_index);
  auto tensor786 = vector<shared_ptr<Tensor>>{I194, t2, I1373};
  auto task786 = make_shared<Task786>(tensor786, pindex);
  task727->add_dep(task786);
  task786->add_dep(task108);
  residualq->add_task(task786);

  auto tensor787 = vector<shared_ptr<Tensor>>{I1373, Gamma48_(), v2_};
  auto task787 = make_shared<Task787>(tensor787, pindex);
  task786->add_dep(task787);
  task787->add_dep(task108);
  residualq->add_task(task787);

  auto tensor788 = vector<shared_ptr<Tensor>>{I1373, Gamma230_(), v2_};
  auto task788 = make_shared<Task788>(tensor788, pindex);
  task786->add_dep(task788);
  task788->add_dep(task108);
  residualq->add_task(task788);

  vector<IndexRange> I1382_index = {virt_, closed_, active_, active_};
  auto I1382 = make_shared<Tensor>(I1382_index);
  auto tensor789 = vector<shared_ptr<Tensor>>{I194, v2_, I1382};
  auto task789 = make_shared<Task789>(tensor789, pindex);
  task727->add_dep(task789);
  task789->add_dep(task108);
  residualq->add_task(task789);

  auto tensor790 = vector<shared_ptr<Tensor>>{I1382, Gamma27_(), t2};
  auto task790 = make_shared<Task790>(tensor790, pindex);
  task789->add_dep(task790);
  task790->add_dep(task108);
  residualq->add_task(task790);

  auto tensor791 = vector<shared_ptr<Tensor>>{I1382, Gamma29_(), t2};
  auto task791 = make_shared<Task791>(tensor791, pindex);
  task789->add_dep(task791);
  task791->add_dep(task108);
  residualq->add_task(task791);

  vector<IndexRange> I1385_index = {virt_, closed_, active_, active_};
  auto I1385 = make_shared<Tensor>(I1385_index);
  auto tensor792 = vector<shared_ptr<Tensor>>{I194, v2_, I1385};
  auto task792 = make_shared<Task792>(tensor792, pindex);
  task727->add_dep(task792);
  task792->add_dep(task108);
  residualq->add_task(task792);

  auto tensor793 = vector<shared_ptr<Tensor>>{I1385, Gamma27_(), t2};
  auto task793 = make_shared<Task793>(tensor793, pindex);
  task792->add_dep(task793);
  task793->add_dep(task108);
  residualq->add_task(task793);

  auto tensor794 = vector<shared_ptr<Tensor>>{I1385, Gamma29_(), t2};
  auto task794 = make_shared<Task794>(tensor794, pindex);
  task792->add_dep(task794);
  task794->add_dep(task108);
  residualq->add_task(task794);

  vector<IndexRange> I1388_index = {virt_, closed_, active_, active_};
  auto I1388 = make_shared<Tensor>(I1388_index);
  auto tensor795 = vector<shared_ptr<Tensor>>{I194, v2_, I1388};
  auto task795 = make_shared<Task795>(tensor795, pindex);
  task727->add_dep(task795);
  task795->add_dep(task108);
  residualq->add_task(task795);

  vector<IndexRange> I1389_index = {active_, virt_, closed_, active_};
  auto I1389 = make_shared<Tensor>(I1389_index);
  auto tensor796 = vector<shared_ptr<Tensor>>{I1388, Gamma29_(), I1389};
  auto task796 = make_shared<Task796>(tensor796, pindex);
  task795->add_dep(task796);
  task796->add_dep(task108);
  residualq->add_task(task796);

  auto tensor797 = vector<shared_ptr<Tensor>>{I1389, t2};
  auto task797 = make_shared<Task797>(tensor797, pindex);
  task796->add_dep(task797);
  task797->add_dep(task108);
  residualq->add_task(task797);

  vector<IndexRange> I1391_index = {virt_, closed_, active_, active_};
  auto I1391 = make_shared<Tensor>(I1391_index);
  auto tensor798 = vector<shared_ptr<Tensor>>{I194, v2_, I1391};
  auto task798 = make_shared<Task798>(tensor798, pindex);
  task727->add_dep(task798);
  task798->add_dep(task108);
  residualq->add_task(task798);

  auto tensor799 = vector<shared_ptr<Tensor>>{I1391, Gamma27_(), t2};
  auto task799 = make_shared<Task799>(tensor799, pindex);
  task798->add_dep(task799);
  task799->add_dep(task108);
  residualq->add_task(task799);

  auto tensor800 = vector<shared_ptr<Tensor>>{I1391, Gamma29_(), t2};
  auto task800 = make_shared<Task800>(tensor800, pindex);
  task798->add_dep(task800);
  task800->add_dep(task108);
  residualq->add_task(task800);

  vector<IndexRange> I1394_index = {virt_, closed_, active_, active_};
  auto I1394 = make_shared<Tensor>(I1394_index);
  auto tensor801 = vector<shared_ptr<Tensor>>{I194, v2_, I1394};
  auto task801 = make_shared<Task801>(tensor801, pindex);
  task727->add_dep(task801);
  task801->add_dep(task108);
  residualq->add_task(task801);

  auto tensor802 = vector<shared_ptr<Tensor>>{I1394, Gamma27_(), t2};
  auto task802 = make_shared<Task802>(tensor802, pindex);
  task801->add_dep(task802);
  task802->add_dep(task108);
  residualq->add_task(task802);

  auto tensor803 = vector<shared_ptr<Tensor>>{I1394, Gamma29_(), t2};
  auto task803 = make_shared<Task803>(tensor803, pindex);
  task801->add_dep(task803);
  task803->add_dep(task108);
  residualq->add_task(task803);

  vector<IndexRange> I1397_index = {virt_, closed_, active_, active_};
  auto I1397 = make_shared<Tensor>(I1397_index);
  auto tensor804 = vector<shared_ptr<Tensor>>{I194, v2_, I1397};
  auto task804 = make_shared<Task804>(tensor804, pindex);
  task727->add_dep(task804);
  task804->add_dep(task108);
  residualq->add_task(task804);

  vector<IndexRange> I1398_index = {active_, virt_, closed_, active_};
  auto I1398 = make_shared<Tensor>(I1398_index);
  auto tensor805 = vector<shared_ptr<Tensor>>{I1397, Gamma29_(), I1398};
  auto task805 = make_shared<Task805>(tensor805, pindex);
  task804->add_dep(task805);
  task805->add_dep(task108);
  residualq->add_task(task805);

  auto tensor806 = vector<shared_ptr<Tensor>>{I1398, t2};
  auto task806 = make_shared<Task806>(tensor806, pindex);
  task805->add_dep(task806);
  task806->add_dep(task108);
  residualq->add_task(task806);

  vector<IndexRange> I1400_index = {virt_, active_, active_, active_};
  auto I1400 = make_shared<Tensor>(I1400_index);
  auto tensor807 = vector<shared_ptr<Tensor>>{I194, t2, I1400};
  auto task807 = make_shared<Task807>(tensor807, pindex);
  task727->add_dep(task807);
  task807->add_dep(task108);
  residualq->add_task(task807);

  auto tensor808 = vector<shared_ptr<Tensor>>{I1400, Gamma49_(), v2_};
  auto task808 = make_shared<Task808>(tensor808, pindex);
  task807->add_dep(task808);
  task808->add_dep(task108);
  residualq->add_task(task808);

  auto tensor809 = vector<shared_ptr<Tensor>>{I1400, Gamma240_(), v2_};
  auto task809 = make_shared<Task809>(tensor809, pindex);
  task807->add_dep(task809);
  task809->add_dep(task108);
  residualq->add_task(task809);

  vector<IndexRange> I1403_index = {virt_, active_, active_, active_};
  auto I1403 = make_shared<Tensor>(I1403_index);
  auto tensor810 = vector<shared_ptr<Tensor>>{I194, t2, I1403};
  auto task810 = make_shared<Task810>(tensor810, pindex);
  task727->add_dep(task810);
  task810->add_dep(task108);
  residualq->add_task(task810);

  auto tensor811 = vector<shared_ptr<Tensor>>{I1403, Gamma49_(), v2_};
  auto task811 = make_shared<Task811>(tensor811, pindex);
  task810->add_dep(task811);
  task811->add_dep(task108);
  residualq->add_task(task811);

  auto tensor812 = vector<shared_ptr<Tensor>>{I1403, Gamma240_(), v2_};
  auto task812 = make_shared<Task812>(tensor812, pindex);
  task810->add_dep(task812);
  task812->add_dep(task108);
  residualq->add_task(task812);

  vector<IndexRange> I1430_index = {virt_, active_, active_, active_};
  auto I1430 = make_shared<Tensor>(I1430_index);
  auto tensor813 = vector<shared_ptr<Tensor>>{I194, v2_, I1430};
  auto task813 = make_shared<Task813>(tensor813, pindex);
  task727->add_dep(task813);
  task813->add_dep(task108);
  residualq->add_task(task813);

  auto tensor814 = vector<shared_ptr<Tensor>>{I1430, Gamma50_(), t2};
  auto task814 = make_shared<Task814>(tensor814, pindex);
  task813->add_dep(task814);
  task814->add_dep(task108);
  residualq->add_task(task814);

  vector<IndexRange> I1433_index = {virt_, active_, active_, active_};
  auto I1433 = make_shared<Tensor>(I1433_index);
  auto tensor815 = vector<shared_ptr<Tensor>>{I194, v2_, I1433};
  auto task815 = make_shared<Task815>(tensor815, pindex);
  task727->add_dep(task815);
  task815->add_dep(task108);
  residualq->add_task(task815);

  auto tensor816 = vector<shared_ptr<Tensor>>{I1433, Gamma50_(), t2};
  auto task816 = make_shared<Task816>(tensor816, pindex);
  task815->add_dep(task816);
  task816->add_dep(task108);
  residualq->add_task(task816);

  vector<IndexRange> I1436_index = {virt_, active_, active_, active_};
  auto I1436 = make_shared<Tensor>(I1436_index);
  auto tensor817 = vector<shared_ptr<Tensor>>{I194, v2_, I1436};
  auto task817 = make_shared<Task817>(tensor817, pindex);
  task727->add_dep(task817);
  task817->add_dep(task108);
  residualq->add_task(task817);

  auto tensor818 = vector<shared_ptr<Tensor>>{I1436, Gamma252_(), t2};
  auto task818 = make_shared<Task818>(tensor818, pindex);
  task817->add_dep(task818);
  task818->add_dep(task108);
  residualq->add_task(task818);

  vector<IndexRange> I1439_index = {virt_, active_, active_, active_};
  auto I1439 = make_shared<Tensor>(I1439_index);
  auto tensor819 = vector<shared_ptr<Tensor>>{I194, v2_, I1439};
  auto task819 = make_shared<Task819>(tensor819, pindex);
  task727->add_dep(task819);
  task819->add_dep(task108);
  residualq->add_task(task819);

  auto tensor820 = vector<shared_ptr<Tensor>>{I1439, Gamma31_(), t2};
  auto task820 = make_shared<Task820>(tensor820, pindex);
  task819->add_dep(task820);
  task820->add_dep(task108);
  residualq->add_task(task820);

  vector<IndexRange> I1442_index = {virt_, active_, active_, active_};
  auto I1442 = make_shared<Tensor>(I1442_index);
  auto tensor821 = vector<shared_ptr<Tensor>>{I194, v2_, I1442};
  auto task821 = make_shared<Task821>(tensor821, pindex);
  task727->add_dep(task821);
  task821->add_dep(task108);
  residualq->add_task(task821);

  auto tensor822 = vector<shared_ptr<Tensor>>{I1442, Gamma471_(), t2};
  auto task822 = make_shared<Task822>(tensor822, pindex);
  task821->add_dep(task822);
  task822->add_dep(task108);
  residualq->add_task(task822);

  vector<IndexRange> I1445_index = {virt_, active_, active_, active_};
  auto I1445 = make_shared<Tensor>(I1445_index);
  auto tensor823 = vector<shared_ptr<Tensor>>{I194, v2_, I1445};
  auto task823 = make_shared<Task823>(tensor823, pindex);
  task727->add_dep(task823);
  task823->add_dep(task108);
  residualq->add_task(task823);

  auto tensor824 = vector<shared_ptr<Tensor>>{I1445, Gamma50_(), t2};
  auto task824 = make_shared<Task824>(tensor824, pindex);
  task823->add_dep(task824);
  task824->add_dep(task108);
  residualq->add_task(task824);

  vector<IndexRange> I1448_index = {virt_, active_, active_, active_};
  auto I1448 = make_shared<Tensor>(I1448_index);
  auto tensor825 = vector<shared_ptr<Tensor>>{I194, v2_, I1448};
  auto task825 = make_shared<Task825>(tensor825, pindex);
  task727->add_dep(task825);
  task825->add_dep(task108);
  residualq->add_task(task825);

  auto tensor826 = vector<shared_ptr<Tensor>>{I1448, Gamma50_(), t2};
  auto task826 = make_shared<Task826>(tensor826, pindex);
  task825->add_dep(task826);
  task826->add_dep(task108);
  residualq->add_task(task826);

  vector<IndexRange> I1451_index = {virt_, active_, active_, active_};
  auto I1451 = make_shared<Tensor>(I1451_index);
  auto tensor827 = vector<shared_ptr<Tensor>>{I194, v2_, I1451};
  auto task827 = make_shared<Task827>(tensor827, pindex);
  task727->add_dep(task827);
  task827->add_dep(task108);
  residualq->add_task(task827);

  auto tensor828 = vector<shared_ptr<Tensor>>{I1451, Gamma50_(), t2};
  auto task828 = make_shared<Task828>(tensor828, pindex);
  task827->add_dep(task828);
  task828->add_dep(task108);
  residualq->add_task(task828);

  vector<IndexRange> I1454_index = {virt_, active_};
  auto I1454 = make_shared<Tensor>(I1454_index);
  auto tensor829 = vector<shared_ptr<Tensor>>{I194, v2_, I1454};
  auto task829 = make_shared<Task829>(tensor829, pindex);
  task727->add_dep(task829);
  task829->add_dep(task108);
  residualq->add_task(task829);

  auto tensor830 = vector<shared_ptr<Tensor>>{I1454, Gamma51_(), t2};
  auto task830 = make_shared<Task830>(tensor830, pindex);
  task829->add_dep(task830);
  task830->add_dep(task108);
  residualq->add_task(task830);

  vector<IndexRange> I1457_index = {virt_, active_};
  auto I1457 = make_shared<Tensor>(I1457_index);
  auto tensor831 = vector<shared_ptr<Tensor>>{I194, v2_, I1457};
  auto task831 = make_shared<Task831>(tensor831, pindex);
  task727->add_dep(task831);
  task831->add_dep(task108);
  residualq->add_task(task831);

  auto tensor832 = vector<shared_ptr<Tensor>>{I1457, Gamma51_(), t2};
  auto task832 = make_shared<Task832>(tensor832, pindex);
  task831->add_dep(task832);
  task832->add_dep(task108);
  residualq->add_task(task832);

  vector<IndexRange> I1472_index = {closed_, closed_, closed_, active_};
  auto I1472 = make_shared<Tensor>(I1472_index);
  auto tensor833 = vector<shared_ptr<Tensor>>{I194, t2, I1472};
  auto task833 = make_shared<Task833>(tensor833, pindex);
  task727->add_dep(task833);
  task833->add_dep(task108);
  residualq->add_task(task833);

  auto tensor834 = vector<shared_ptr<Tensor>>{I1472, Gamma32_(), v2_};
  auto task834 = make_shared<Task834>(tensor834, pindex);
  task833->add_dep(task834);
  task834->add_dep(task108);
  residualq->add_task(task834);

  vector<IndexRange> I1475_index = {closed_, closed_, closed_, active_};
  auto I1475 = make_shared<Tensor>(I1475_index);
  auto tensor835 = vector<shared_ptr<Tensor>>{I194, t2, I1475};
  auto task835 = make_shared<Task835>(tensor835, pindex);
  task727->add_dep(task835);
  task835->add_dep(task108);
  residualq->add_task(task835);

  auto tensor836 = vector<shared_ptr<Tensor>>{I1475, Gamma32_(), v2_};
  auto task836 = make_shared<Task836>(tensor836, pindex);
  task835->add_dep(task836);
  task836->add_dep(task108);
  residualq->add_task(task836);

  vector<IndexRange> I1502_index = {closed_, closed_, active_, active_};
  auto I1502 = make_shared<Tensor>(I1502_index);
  auto tensor837 = vector<shared_ptr<Tensor>>{I194, t2, I1502};
  auto task837 = make_shared<Task837>(tensor837, pindex);
  task727->add_dep(task837);
  task837->add_dep(task108);
  residualq->add_task(task837);

  vector<IndexRange> I1503_index = {closed_, closed_, active_, active_};
  auto I1503 = make_shared<Tensor>(I1503_index);
  auto tensor838 = vector<shared_ptr<Tensor>>{I1502, Gamma51_(), I1503};
  auto task838 = make_shared<Task838>(tensor838, pindex);
  task837->add_dep(task838);
  task838->add_dep(task108);
  residualq->add_task(task838);

  auto tensor839 = vector<shared_ptr<Tensor>>{I1503, v2_};
  auto task839 = make_shared<Task839>(tensor839, pindex);
  task838->add_dep(task839);
  task839->add_dep(task108);
  residualq->add_task(task839);

  auto tensor840 = vector<shared_ptr<Tensor>>{I1502, Gamma29_(), v2_};
  auto task840 = make_shared<Task840>(tensor840, pindex);
  task837->add_dep(task840);
  task840->add_dep(task108);
  residualq->add_task(task840);

  auto tensor841 = vector<shared_ptr<Tensor>>{I1502, Gamma503_(), v2_};
  auto task841 = make_shared<Task841>(tensor841, pindex);
  task837->add_dep(task841);
  task841->add_dep(task108);
  residualq->add_task(task841);

  vector<IndexRange> I1505_index = {closed_, closed_, active_, active_};
  auto I1505 = make_shared<Tensor>(I1505_index);
  auto tensor842 = vector<shared_ptr<Tensor>>{I194, t2, I1505};
  auto task842 = make_shared<Task842>(tensor842, pindex);
  task727->add_dep(task842);
  task842->add_dep(task108);
  residualq->add_task(task842);

  vector<IndexRange> I1506_index = {closed_, closed_, active_, active_};
  auto I1506 = make_shared<Tensor>(I1506_index);
  auto tensor843 = vector<shared_ptr<Tensor>>{I1505, Gamma51_(), I1506};
  auto task843 = make_shared<Task843>(tensor843, pindex);
  task842->add_dep(task843);
  task843->add_dep(task108);
  residualq->add_task(task843);

  auto tensor844 = vector<shared_ptr<Tensor>>{I1506, v2_};
  auto task844 = make_shared<Task844>(tensor844, pindex);
  task843->add_dep(task844);
  task844->add_dep(task108);
  residualq->add_task(task844);

  auto tensor845 = vector<shared_ptr<Tensor>>{I1505, Gamma27_(), v2_};
  auto task845 = make_shared<Task845>(tensor845, pindex);
  task842->add_dep(task845);
  task845->add_dep(task108);
  residualq->add_task(task845);

  vector<IndexRange> I1508_index = {virt_, virt_, active_, active_};
  auto I1508 = make_shared<Tensor>(I1508_index);
  auto tensor846 = vector<shared_ptr<Tensor>>{I194, t2, I1508};
  auto task846 = make_shared<Task846>(tensor846, pindex);
  task727->add_dep(task846);
  task846->add_dep(task108);
  residualq->add_task(task846);

  vector<IndexRange> I1509_index = {virt_, virt_, active_, active_};
  auto I1509 = make_shared<Tensor>(I1509_index);
  auto tensor847 = vector<shared_ptr<Tensor>>{I1508, Gamma51_(), I1509};
  auto task847 = make_shared<Task847>(tensor847, pindex);
  task846->add_dep(task847);
  task847->add_dep(task108);
  residualq->add_task(task847);

  auto tensor848 = vector<shared_ptr<Tensor>>{I1509, v2_};
  auto task848 = make_shared<Task848>(tensor848, pindex);
  task847->add_dep(task848);
  task848->add_dep(task108);
  residualq->add_task(task848);

  auto tensor849 = vector<shared_ptr<Tensor>>{I1508, Gamma29_(), v2_};
  auto task849 = make_shared<Task849>(tensor849, pindex);
  task846->add_dep(task849);
  task849->add_dep(task108);
  residualq->add_task(task849);

  auto tensor850 = vector<shared_ptr<Tensor>>{I1508, Gamma503_(), v2_};
  auto task850 = make_shared<Task850>(tensor850, pindex);
  task846->add_dep(task850);
  task850->add_dep(task108);
  residualq->add_task(task850);

  vector<IndexRange> I1511_index = {virt_, virt_, active_, active_};
  auto I1511 = make_shared<Tensor>(I1511_index);
  auto tensor851 = vector<shared_ptr<Tensor>>{I194, t2, I1511};
  auto task851 = make_shared<Task851>(tensor851, pindex);
  task727->add_dep(task851);
  task851->add_dep(task108);
  residualq->add_task(task851);

  vector<IndexRange> I1512_index = {virt_, virt_, active_, active_};
  auto I1512 = make_shared<Tensor>(I1512_index);
  auto tensor852 = vector<shared_ptr<Tensor>>{I1511, Gamma51_(), I1512};
  auto task852 = make_shared<Task852>(tensor852, pindex);
  task851->add_dep(task852);
  task852->add_dep(task108);
  residualq->add_task(task852);

  auto tensor853 = vector<shared_ptr<Tensor>>{I1512, v2_};
  auto task853 = make_shared<Task853>(tensor853, pindex);
  task852->add_dep(task853);
  task853->add_dep(task108);
  residualq->add_task(task853);

  auto tensor854 = vector<shared_ptr<Tensor>>{I1511, Gamma29_(), v2_};
  auto task854 = make_shared<Task854>(tensor854, pindex);
  task851->add_dep(task854);
  task854->add_dep(task108);
  residualq->add_task(task854);

  auto tensor855 = vector<shared_ptr<Tensor>>{I1511, Gamma503_(), v2_};
  auto task855 = make_shared<Task855>(tensor855, pindex);
  task851->add_dep(task855);
  task855->add_dep(task108);
  residualq->add_task(task855);

  vector<IndexRange> I1514_index = {virt_, virt_, active_, active_};
  auto I1514 = make_shared<Tensor>(I1514_index);
  auto tensor856 = vector<shared_ptr<Tensor>>{I194, t2, I1514};
  auto task856 = make_shared<Task856>(tensor856, pindex);
  task727->add_dep(task856);
  task856->add_dep(task108);
  residualq->add_task(task856);

  vector<IndexRange> I1515_index = {virt_, virt_, active_, active_};
  auto I1515 = make_shared<Tensor>(I1515_index);
  auto tensor857 = vector<shared_ptr<Tensor>>{I1514, Gamma51_(), I1515};
  auto task857 = make_shared<Task857>(tensor857, pindex);
  task856->add_dep(task857);
  task857->add_dep(task108);
  residualq->add_task(task857);

  auto tensor858 = vector<shared_ptr<Tensor>>{I1515, v2_};
  auto task858 = make_shared<Task858>(tensor858, pindex);
  task857->add_dep(task858);
  task858->add_dep(task108);
  residualq->add_task(task858);

  auto tensor859 = vector<shared_ptr<Tensor>>{I1514, Gamma29_(), v2_};
  auto task859 = make_shared<Task859>(tensor859, pindex);
  task856->add_dep(task859);
  task859->add_dep(task108);
  residualq->add_task(task859);

  auto tensor860 = vector<shared_ptr<Tensor>>{I1514, Gamma503_(), v2_};
  auto task860 = make_shared<Task860>(tensor860, pindex);
  task856->add_dep(task860);
  task860->add_dep(task108);
  residualq->add_task(task860);

  vector<IndexRange> I1517_index = {virt_, virt_, active_, active_};
  auto I1517 = make_shared<Tensor>(I1517_index);
  auto tensor861 = vector<shared_ptr<Tensor>>{I194, t2, I1517};
  auto task861 = make_shared<Task861>(tensor861, pindex);
  task727->add_dep(task861);
  task861->add_dep(task108);
  residualq->add_task(task861);

  vector<IndexRange> I1518_index = {virt_, virt_, active_, active_};
  auto I1518 = make_shared<Tensor>(I1518_index);
  auto tensor862 = vector<shared_ptr<Tensor>>{I1517, Gamma51_(), I1518};
  auto task862 = make_shared<Task862>(tensor862, pindex);
  task861->add_dep(task862);
  task862->add_dep(task108);
  residualq->add_task(task862);

  auto tensor863 = vector<shared_ptr<Tensor>>{I1518, v2_};
  auto task863 = make_shared<Task863>(tensor863, pindex);
  task862->add_dep(task863);
  task863->add_dep(task108);
  residualq->add_task(task863);

  auto tensor864 = vector<shared_ptr<Tensor>>{I1517, Gamma27_(), v2_};
  auto task864 = make_shared<Task864>(tensor864, pindex);
  task861->add_dep(task864);
  task864->add_dep(task108);
  residualq->add_task(task864);

  vector<IndexRange> I1604_index = {closed_, active_, active_, active_};
  auto I1604 = make_shared<Tensor>(I1604_index);
  auto tensor865 = vector<shared_ptr<Tensor>>{I194, t2, I1604};
  auto task865 = make_shared<Task865>(tensor865, pindex);
  task727->add_dep(task865);
  task865->add_dep(task108);
  residualq->add_task(task865);

  auto tensor866 = vector<shared_ptr<Tensor>>{I1604, Gamma50_(), v2_};
  auto task866 = make_shared<Task866>(tensor866, pindex);
  task865->add_dep(task866);
  task866->add_dep(task108);
  residualq->add_task(task866);

  auto tensor867 = vector<shared_ptr<Tensor>>{I1604, Gamma526_(), v2_};
  auto task867 = make_shared<Task867>(tensor867, pindex);
  task865->add_dep(task867);
  task867->add_dep(task108);
  residualq->add_task(task867);

  vector<IndexRange> I1610_index = {virt_, virt_, active_, active_};
  auto I1610 = make_shared<Tensor>(I1610_index);
  auto tensor868 = vector<shared_ptr<Tensor>>{I194, v2_, I1610};
  auto task868 = make_shared<Task868>(tensor868, pindex);
  task727->add_dep(task868);
  task868->add_dep(task108);
  residualq->add_task(task868);

  auto tensor869 = vector<shared_ptr<Tensor>>{I1610, Gamma51_(), t2};
  auto task869 = make_shared<Task869>(tensor869, pindex);
  task868->add_dep(task869);
  task869->add_dep(task108);
  residualq->add_task(task869);

  vector<IndexRange> I1613_index = {virt_, virt_, active_, active_};
  auto I1613 = make_shared<Tensor>(I1613_index);
  auto tensor870 = vector<shared_ptr<Tensor>>{I194, v2_, I1613};
  auto task870 = make_shared<Task870>(tensor870, pindex);
  task727->add_dep(task870);
  task870->add_dep(task108);
  residualq->add_task(task870);

  auto tensor871 = vector<shared_ptr<Tensor>>{I1613, Gamma51_(), t2};
  auto task871 = make_shared<Task871>(tensor871, pindex);
  task870->add_dep(task871);
  task871->add_dep(task108);
  residualq->add_task(task871);

  vector<IndexRange> I1616_index = {virt_, virt_, active_, active_};
  auto I1616 = make_shared<Tensor>(I1616_index);
  auto tensor872 = vector<shared_ptr<Tensor>>{I194, v2_, I1616};
  auto task872 = make_shared<Task872>(tensor872, pindex);
  task727->add_dep(task872);
  task872->add_dep(task108);
  residualq->add_task(task872);

  auto tensor873 = vector<shared_ptr<Tensor>>{I1616, Gamma503_(), t2};
  auto task873 = make_shared<Task873>(tensor873, pindex);
  task872->add_dep(task873);
  task873->add_dep(task108);
  residualq->add_task(task873);

  vector<IndexRange> I1619_index = {virt_, virt_, active_, active_};
  auto I1619 = make_shared<Tensor>(I1619_index);
  auto tensor874 = vector<shared_ptr<Tensor>>{I194, v2_, I1619};
  auto task874 = make_shared<Task874>(tensor874, pindex);
  task727->add_dep(task874);
  task874->add_dep(task108);
  residualq->add_task(task874);

  auto tensor875 = vector<shared_ptr<Tensor>>{I1619, Gamma51_(), t2};
  auto task875 = make_shared<Task875>(tensor875, pindex);
  task874->add_dep(task875);
  task875->add_dep(task108);
  residualq->add_task(task875);

  vector<IndexRange> I1701_index = {active_, virt_, closed_, virt_};
  auto I1701 = make_shared<Tensor>(I1701_index);
  auto tensor876 = vector<shared_ptr<Tensor>>{I194, Gamma562_(), I1701};
  auto task876 = make_shared<Task876>(tensor876, pindex);
  task727->add_dep(task876);
  task876->add_dep(task108);
  residualq->add_task(task876);

  auto tensor877 = vector<shared_ptr<Tensor>>{I1701, t2};
  auto task877 = make_shared<Task877>(tensor877, pindex);
  task876->add_dep(task877);
  task877->add_dep(task108);
  residualq->add_task(task877);

  vector<IndexRange> I1705_index = {active_, virt_, closed_, virt_};
  auto I1705 = make_shared<Tensor>(I1705_index);
  auto tensor878 = vector<shared_ptr<Tensor>>{I194, Gamma564_(), I1705};
  auto task878 = make_shared<Task878>(tensor878, pindex);
  task727->add_dep(task878);
  task878->add_dep(task108);
  residualq->add_task(task878);

  auto tensor879 = vector<shared_ptr<Tensor>>{I1705, t2};
  auto task879 = make_shared<Task879>(tensor879, pindex);
  task878->add_dep(task879);
  task879->add_dep(task108);
  residualq->add_task(task879);
}

#endif
