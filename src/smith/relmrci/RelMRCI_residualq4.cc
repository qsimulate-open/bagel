//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq4.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks9.h>
#include <src/smith/relmrci/RelMRCI_tasks10.h>
#include <src/smith/relmrci/RelMRCI_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void RelMRCI::RelMRCI::make_residualq4(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I108_index = {virt_, closed_, virt_, closed_};
  auto I108 = make_shared<Tensor>(I108_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{r, I108};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task426->add_dep(task83);
  residualq->add_task(task426);

  vector<IndexRange> I109_index = {virt_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor427 = vector<shared_ptr<Tensor>>{I108, t2, I109};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task83);
  residualq->add_task(task427);

  auto tensor428 = vector<shared_ptr<Tensor>>{I109, Gamma11_(), h1_};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task427->add_dep(task428);
  task428->add_dep(task83);
  residualq->add_task(task428);

  vector<IndexRange> I112_index = {virt_, active_};
  auto I112 = make_shared<Tensor>(I112_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{I108, t2, I112};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task426->add_dep(task429);
  task429->add_dep(task83);
  residualq->add_task(task429);

  auto tensor430 = vector<shared_ptr<Tensor>>{I112, Gamma11_(), h1_};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task83);
  residualq->add_task(task430);

  vector<IndexRange> I115_index = {closed_, virt_};
  auto I115 = make_shared<Tensor>(I115_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{I108, h1_, I115};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task426->add_dep(task431);
  task431->add_dep(task83);
  residualq->add_task(task431);

  auto tensor432 = vector<shared_ptr<Tensor>>{I115, Gamma27_(), t2};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task83);
  residualq->add_task(task432);

  vector<IndexRange> I118_index = {closed_, virt_};
  auto I118 = make_shared<Tensor>(I118_index);
  auto tensor433 = vector<shared_ptr<Tensor>>{I108, h1_, I118};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task426->add_dep(task433);
  task433->add_dep(task83);
  residualq->add_task(task433);

  auto tensor434 = vector<shared_ptr<Tensor>>{I118, Gamma27_(), t2};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task433->add_dep(task434);
  task434->add_dep(task83);
  residualq->add_task(task434);

  vector<IndexRange> I121_index = {closed_, closed_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I108, t2, I121};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task426->add_dep(task435);
  task435->add_dep(task83);
  residualq->add_task(task435);

  shared_ptr<Task436> task436;
  if (diagonal) {
    auto tensor436 = vector<shared_ptr<Tensor>>{I121, h1_};
    task436 = make_shared<Task436>(tensor436, pindex);
    task435->add_dep(task436);
    task436->add_dep(task83);
    residualq->add_task(task436);
  }

  vector<IndexRange> I868_index = {closed_, closed_, active_, active_};
  auto I868 = make_shared<Tensor>(I868_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{I121, Gamma27_(), I868};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task435->add_dep(task437);
  task437->add_dep(task83);
  residualq->add_task(task437);

  auto tensor438 = vector<shared_ptr<Tensor>>{I868, v2_};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task83);
  residualq->add_task(task438);

  auto tensor439 = vector<shared_ptr<Tensor>>{I121, Gamma11_(), v2_};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task435->add_dep(task439);
  task439->add_dep(task83);
  residualq->add_task(task439);

  vector<IndexRange> I123_index = {closed_, closed_};
  auto I123 = make_shared<Tensor>(I123_index);
  auto tensor440 = vector<shared_ptr<Tensor>>{I108, t2, I123};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task426->add_dep(task440);
  task440->add_dep(task83);
  residualq->add_task(task440);

  shared_ptr<Task441> task441;
  if (diagonal) {
    auto tensor441 = vector<shared_ptr<Tensor>>{I123, h1_};
    task441 = make_shared<Task441>(tensor441, pindex);
    task440->add_dep(task441);
    task441->add_dep(task83);
    residualq->add_task(task441);
  }

  vector<IndexRange> I871_index = {closed_, closed_, active_, active_};
  auto I871 = make_shared<Tensor>(I871_index);
  auto tensor442 = vector<shared_ptr<Tensor>>{I123, Gamma27_(), I871};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task440->add_dep(task442);
  task442->add_dep(task83);
  residualq->add_task(task442);

  auto tensor443 = vector<shared_ptr<Tensor>>{I871, v2_};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task442->add_dep(task443);
  task443->add_dep(task83);
  residualq->add_task(task443);

  auto tensor444 = vector<shared_ptr<Tensor>>{I123, Gamma11_(), v2_};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task440->add_dep(task444);
  task444->add_dep(task83);
  residualq->add_task(task444);

  vector<IndexRange> I125_index = {virt_, virt_};
  auto I125 = make_shared<Tensor>(I125_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I108, t2, I125};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task426->add_dep(task445);
  task445->add_dep(task83);
  residualq->add_task(task445);

  shared_ptr<Task446> task446;
  if (diagonal) {
    auto tensor446 = vector<shared_ptr<Tensor>>{I125, h1_};
    task446 = make_shared<Task446>(tensor446, pindex);
    task445->add_dep(task446);
    task446->add_dep(task83);
    residualq->add_task(task446);
  }

  vector<IndexRange> I874_index = {virt_, virt_, active_, active_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{I125, Gamma27_(), I874};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task445->add_dep(task447);
  task447->add_dep(task83);
  residualq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I874, v2_};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task83);
  residualq->add_task(task448);

  auto tensor449 = vector<shared_ptr<Tensor>>{I125, Gamma11_(), v2_};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task445->add_dep(task449);
  task449->add_dep(task83);
  residualq->add_task(task449);

  shared_ptr<Task450> task450;
  if (diagonal) {
    auto tensor450 = vector<shared_ptr<Tensor>>{I108, h1_, t2};
    task450 = make_shared<Task450>(tensor450, pindex);
    task426->add_dep(task450);
    task450->add_dep(task83);
    residualq->add_task(task450);
  }

  vector<IndexRange> I129_index = {closed_, active_};
  auto I129 = make_shared<Tensor>(I129_index);
  auto tensor451 = vector<shared_ptr<Tensor>>{I108, t2, I129};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task426->add_dep(task451);
  task451->add_dep(task83);
  residualq->add_task(task451);

  auto tensor452 = vector<shared_ptr<Tensor>>{I129, Gamma27_(), h1_};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task451->add_dep(task452);
  task452->add_dep(task83);
  residualq->add_task(task452);

  vector<IndexRange> I132_index = {closed_, active_};
  auto I132 = make_shared<Tensor>(I132_index);
  auto tensor453 = vector<shared_ptr<Tensor>>{I108, t2, I132};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task426->add_dep(task453);
  task453->add_dep(task83);
  residualq->add_task(task453);

  auto tensor454 = vector<shared_ptr<Tensor>>{I132, Gamma27_(), h1_};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task453->add_dep(task454);
  task454->add_dep(task83);
  residualq->add_task(task454);

  vector<IndexRange> I777_index = {closed_, active_};
  auto I777 = make_shared<Tensor>(I777_index);
  auto tensor455 = vector<shared_ptr<Tensor>>{I108, v2_, I777};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task426->add_dep(task455);
  task455->add_dep(task83);
  residualq->add_task(task455);

  auto tensor456 = vector<shared_ptr<Tensor>>{I777, Gamma9_(), t2};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task455->add_dep(task456);
  task456->add_dep(task83);
  residualq->add_task(task456);

  vector<IndexRange> I780_index = {closed_, active_};
  auto I780 = make_shared<Tensor>(I780_index);
  auto tensor457 = vector<shared_ptr<Tensor>>{I108, v2_, I780};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task426->add_dep(task457);
  task457->add_dep(task83);
  residualq->add_task(task457);

  auto tensor458 = vector<shared_ptr<Tensor>>{I780, Gamma9_(), t2};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task83);
  residualq->add_task(task458);

  vector<IndexRange> I783_index = {virt_, active_};
  auto I783 = make_shared<Tensor>(I783_index);
  auto tensor459 = vector<shared_ptr<Tensor>>{I108, t2, I783};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task426->add_dep(task459);
  task459->add_dep(task83);
  residualq->add_task(task459);

  auto tensor460 = vector<shared_ptr<Tensor>>{I783, Gamma5_(), v2_};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task459->add_dep(task460);
  task460->add_dep(task83);
  residualq->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I783, Gamma160_(), v2_};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task459->add_dep(task461);
  task461->add_dep(task83);
  residualq->add_task(task461);

  vector<IndexRange> I786_index = {virt_, active_};
  auto I786 = make_shared<Tensor>(I786_index);
  auto tensor462 = vector<shared_ptr<Tensor>>{I108, t2, I786};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task426->add_dep(task462);
  task462->add_dep(task83);
  residualq->add_task(task462);

  auto tensor463 = vector<shared_ptr<Tensor>>{I786, Gamma5_(), v2_};
  auto task463 = make_shared<Task463>(tensor463, pindex);
  task462->add_dep(task463);
  task463->add_dep(task83);
  residualq->add_task(task463);

  auto tensor464 = vector<shared_ptr<Tensor>>{I786, Gamma160_(), v2_};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task462->add_dep(task464);
  task464->add_dep(task83);
  residualq->add_task(task464);

  vector<IndexRange> I795_index = {closed_, closed_, virt_, active_};
  auto I795 = make_shared<Tensor>(I795_index);
  auto tensor465 = vector<shared_ptr<Tensor>>{I108, t2, I795};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task426->add_dep(task465);
  task465->add_dep(task83);
  residualq->add_task(task465);

  vector<IndexRange> I796_index = {active_, closed_, closed_, virt_};
  auto I796 = make_shared<Tensor>(I796_index);
  auto tensor466 = vector<shared_ptr<Tensor>>{I795, Gamma11_(), I796};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task465->add_dep(task466);
  task466->add_dep(task83);
  residualq->add_task(task466);

  auto tensor467 = vector<shared_ptr<Tensor>>{I796, v2_};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task466->add_dep(task467);
  task467->add_dep(task83);
  residualq->add_task(task467);

  vector<IndexRange> I798_index = {closed_, closed_, virt_, active_};
  auto I798 = make_shared<Tensor>(I798_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{I108, t2, I798};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task426->add_dep(task468);
  task468->add_dep(task83);
  residualq->add_task(task468);

  vector<IndexRange> I799_index = {active_, closed_, closed_, virt_};
  auto I799 = make_shared<Tensor>(I799_index);
  auto tensor469 = vector<shared_ptr<Tensor>>{I798, Gamma11_(), I799};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task83);
  residualq->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I799, v2_};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task469->add_dep(task470);
  task470->add_dep(task83);
  residualq->add_task(task470);

  vector<IndexRange> I807_index = {closed_, closed_, virt_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor471 = vector<shared_ptr<Tensor>>{I108, t2, I807};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task426->add_dep(task471);
  task471->add_dep(task83);
  residualq->add_task(task471);

  vector<IndexRange> I808_index = {active_, closed_, closed_, virt_};
  auto I808 = make_shared<Tensor>(I808_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{I807, Gamma11_(), I808};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task471->add_dep(task472);
  task472->add_dep(task83);
  residualq->add_task(task472);

  auto tensor473 = vector<shared_ptr<Tensor>>{I808, v2_};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task472->add_dep(task473);
  task473->add_dep(task83);
  residualq->add_task(task473);

  vector<IndexRange> I810_index = {closed_, closed_, virt_, active_};
  auto I810 = make_shared<Tensor>(I810_index);
  auto tensor474 = vector<shared_ptr<Tensor>>{I108, t2, I810};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task426->add_dep(task474);
  task474->add_dep(task83);
  residualq->add_task(task474);

  vector<IndexRange> I811_index = {active_, closed_, closed_, virt_};
  auto I811 = make_shared<Tensor>(I811_index);
  auto tensor475 = vector<shared_ptr<Tensor>>{I810, Gamma11_(), I811};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task474->add_dep(task475);
  task475->add_dep(task83);
  residualq->add_task(task475);

  auto tensor476 = vector<shared_ptr<Tensor>>{I811, v2_};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task475->add_dep(task476);
  task476->add_dep(task83);
  residualq->add_task(task476);

  vector<IndexRange> I819_index = {closed_, virt_, closed_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor477 = vector<shared_ptr<Tensor>>{I108, v2_, I819};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task426->add_dep(task477);
  task477->add_dep(task83);
  residualq->add_task(task477);

  auto tensor478 = vector<shared_ptr<Tensor>>{I819, Gamma11_(), t2};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task477->add_dep(task478);
  task478->add_dep(task83);
  residualq->add_task(task478);

  vector<IndexRange> I822_index = {closed_, virt_, closed_, active_};
  auto I822 = make_shared<Tensor>(I822_index);
  auto tensor479 = vector<shared_ptr<Tensor>>{I108, v2_, I822};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task426->add_dep(task479);
  task479->add_dep(task83);
  residualq->add_task(task479);

  auto tensor480 = vector<shared_ptr<Tensor>>{I822, Gamma11_(), t2};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task479->add_dep(task480);
  task480->add_dep(task83);
  residualq->add_task(task480);

  vector<IndexRange> I825_index = {closed_, virt_, active_, active_};
  auto I825 = make_shared<Tensor>(I825_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{I108, t2, I825};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task426->add_dep(task481);
  task481->add_dep(task83);
  residualq->add_task(task481);

  vector<IndexRange> I826_index = {closed_, virt_, active_, active_};
  auto I826 = make_shared<Tensor>(I826_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{I825, Gamma24_(), I826};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task481->add_dep(task482);
  task482->add_dep(task83);
  residualq->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I826, v2_};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task83);
  residualq->add_task(task483);

  auto tensor484 = vector<shared_ptr<Tensor>>{I825, Gamma9_(), v2_};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task481->add_dep(task484);
  task484->add_dep(task83);
  residualq->add_task(task484);

  vector<IndexRange> I828_index = {closed_, virt_, active_, active_};
  auto I828 = make_shared<Tensor>(I828_index);
  auto tensor485 = vector<shared_ptr<Tensor>>{I108, t2, I828};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task426->add_dep(task485);
  task485->add_dep(task83);
  residualq->add_task(task485);

  vector<IndexRange> I829_index = {closed_, virt_, active_, active_};
  auto I829 = make_shared<Tensor>(I829_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{I828, Gamma24_(), I829};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task485->add_dep(task486);
  task486->add_dep(task83);
  residualq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I829, v2_};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task83);
  residualq->add_task(task487);

  auto tensor488 = vector<shared_ptr<Tensor>>{I828, Gamma9_(), v2_};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task485->add_dep(task488);
  task488->add_dep(task83);
  residualq->add_task(task488);

  vector<IndexRange> I849_index = {closed_, virt_};
  auto I849 = make_shared<Tensor>(I849_index);
  auto tensor489 = vector<shared_ptr<Tensor>>{I108, v2_, I849};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task426->add_dep(task489);
  task489->add_dep(task83);
  residualq->add_task(task489);

  auto tensor490 = vector<shared_ptr<Tensor>>{I849, Gamma27_(), t2};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task489->add_dep(task490);
  task490->add_dep(task83);
  residualq->add_task(task490);

  vector<IndexRange> I852_index = {closed_, virt_};
  auto I852 = make_shared<Tensor>(I852_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{I108, v2_, I852};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task426->add_dep(task491);
  task491->add_dep(task83);
  residualq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I852, Gamma27_(), t2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task83);
  residualq->add_task(task492);

  vector<IndexRange> I855_index = {closed_, virt_};
  auto I855 = make_shared<Tensor>(I855_index);
  auto tensor493 = vector<shared_ptr<Tensor>>{I108, v2_, I855};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task426->add_dep(task493);
  task493->add_dep(task83);
  residualq->add_task(task493);

  auto tensor494 = vector<shared_ptr<Tensor>>{I855, Gamma27_(), t2};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task493->add_dep(task494);
  task494->add_dep(task83);
  residualq->add_task(task494);

  vector<IndexRange> I858_index = {closed_, virt_};
  auto I858 = make_shared<Tensor>(I858_index);
  auto tensor495 = vector<shared_ptr<Tensor>>{I108, v2_, I858};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task426->add_dep(task495);
  task495->add_dep(task83);
  residualq->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I858, Gamma27_(), t2};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task495->add_dep(task496);
  task496->add_dep(task83);
  residualq->add_task(task496);

  vector<IndexRange> I861_index = {virt_, active_};
  auto I861 = make_shared<Tensor>(I861_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{I108, v2_, I861};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task426->add_dep(task497);
  task497->add_dep(task83);
  residualq->add_task(task497);

  auto tensor498 = vector<shared_ptr<Tensor>>{I861, Gamma33_(), t2};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task83);
  residualq->add_task(task498);

  vector<IndexRange> I864_index = {virt_, active_};
  auto I864 = make_shared<Tensor>(I864_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{I108, v2_, I864};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task426->add_dep(task499);
  task499->add_dep(task83);
  residualq->add_task(task499);

  auto tensor500 = vector<shared_ptr<Tensor>>{I864, Gamma33_(), t2};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task83);
  residualq->add_task(task500);

  vector<IndexRange> I876_index = {virt_, virt_};
  auto I876 = make_shared<Tensor>(I876_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{I108, t2, I876};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task426->add_dep(task501);
  task501->add_dep(task83);
  residualq->add_task(task501);

  vector<IndexRange> I877_index = {virt_, virt_, active_, active_};
  auto I877 = make_shared<Tensor>(I877_index);
  auto tensor502 = vector<shared_ptr<Tensor>>{I876, Gamma27_(), I877};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task83);
  residualq->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I877, v2_};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task502->add_dep(task503);
  task503->add_dep(task83);
  residualq->add_task(task503);

  auto tensor504 = vector<shared_ptr<Tensor>>{I876, Gamma11_(), v2_};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task501->add_dep(task504);
  task504->add_dep(task83);
  residualq->add_task(task504);

  shared_ptr<Task505> task505;
  if (diagonal) {
    auto tensor505 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task505 = make_shared<Task505>(tensor505, pindex);
    task426->add_dep(task505);
    task505->add_dep(task83);
    residualq->add_task(task505);
  }

  shared_ptr<Task506> task506;
  if (diagonal) {
    auto tensor506 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task506 = make_shared<Task506>(tensor506, pindex);
    task426->add_dep(task506);
    task506->add_dep(task83);
    residualq->add_task(task506);
  }

  shared_ptr<Task507> task507;
  if (diagonal) {
    auto tensor507 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task507 = make_shared<Task507>(tensor507, pindex);
    task426->add_dep(task507);
    task507->add_dep(task83);
    residualq->add_task(task507);
  }

  shared_ptr<Task508> task508;
  if (diagonal) {
    auto tensor508 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task508 = make_shared<Task508>(tensor508, pindex);
    task426->add_dep(task508);
    task508->add_dep(task83);
    residualq->add_task(task508);
  }

  shared_ptr<Task509> task509;
  if (diagonal) {
    auto tensor509 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task509 = make_shared<Task509>(tensor509, pindex);
    task426->add_dep(task509);
    task509->add_dep(task83);
    residualq->add_task(task509);
  }

  shared_ptr<Task510> task510;
  if (diagonal) {
    auto tensor510 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task510 = make_shared<Task510>(tensor510, pindex);
    task426->add_dep(task510);
    task510->add_dep(task83);
    residualq->add_task(task510);
  }

  shared_ptr<Task511> task511;
  if (diagonal) {
    auto tensor511 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task511 = make_shared<Task511>(tensor511, pindex);
    task426->add_dep(task511);
    task511->add_dep(task83);
    residualq->add_task(task511);
  }

  shared_ptr<Task512> task512;
  if (diagonal) {
    auto tensor512 = vector<shared_ptr<Tensor>>{I108, t2, v2_};
    task512 = make_shared<Task512>(tensor512, pindex);
    task426->add_dep(task512);
    task512->add_dep(task83);
    residualq->add_task(task512);
  }

  vector<IndexRange> I935_index = {closed_, active_};
  auto I935 = make_shared<Tensor>(I935_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{I108, t2, I935};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task426->add_dep(task513);
  task513->add_dep(task83);
  residualq->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I935, Gamma24_(), v2_};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task83);
  residualq->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I935, Gamma33_(), v2_};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task513->add_dep(task515);
  task515->add_dep(task83);
  residualq->add_task(task515);

  vector<IndexRange> I938_index = {closed_, active_};
  auto I938 = make_shared<Tensor>(I938_index);
  auto tensor516 = vector<shared_ptr<Tensor>>{I108, t2, I938};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task426->add_dep(task516);
  task516->add_dep(task83);
  residualq->add_task(task516);

  auto tensor517 = vector<shared_ptr<Tensor>>{I938, Gamma24_(), v2_};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task516->add_dep(task517);
  task517->add_dep(task83);
  residualq->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I938, Gamma33_(), v2_};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task516->add_dep(task518);
  task518->add_dep(task83);
  residualq->add_task(task518);

  vector<IndexRange> I947_index = {closed_, closed_, closed_, active_};
  auto I947 = make_shared<Tensor>(I947_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{I108, t2, I947};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task426->add_dep(task519);
  task519->add_dep(task83);
  residualq->add_task(task519);

  auto tensor520 = vector<shared_ptr<Tensor>>{I947, Gamma27_(), v2_};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task519->add_dep(task520);
  task520->add_dep(task83);
  residualq->add_task(task520);

  vector<IndexRange> I950_index = {closed_, closed_, closed_, active_};
  auto I950 = make_shared<Tensor>(I950_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{I108, t2, I950};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task426->add_dep(task521);
  task521->add_dep(task83);
  residualq->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I950, Gamma27_(), v2_};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task83);
  residualq->add_task(task522);

  vector<IndexRange> I953_index = {virt_, closed_, virt_, active_};
  auto I953 = make_shared<Tensor>(I953_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{I108, t2, I953};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task426->add_dep(task523);
  task523->add_dep(task83);
  residualq->add_task(task523);

  vector<IndexRange> I954_index = {virt_, active_, closed_, virt_};
  auto I954 = make_shared<Tensor>(I954_index);
  auto tensor524 = vector<shared_ptr<Tensor>>{I953, Gamma27_(), I954};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task83);
  residualq->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I954, v2_};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task83);
  residualq->add_task(task525);

  vector<IndexRange> I956_index = {virt_, closed_, virt_, active_};
  auto I956 = make_shared<Tensor>(I956_index);
  auto tensor526 = vector<shared_ptr<Tensor>>{I108, t2, I956};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task426->add_dep(task526);
  task526->add_dep(task83);
  residualq->add_task(task526);

  vector<IndexRange> I957_index = {virt_, active_, closed_, virt_};
  auto I957 = make_shared<Tensor>(I957_index);
  auto tensor527 = vector<shared_ptr<Tensor>>{I956, Gamma27_(), I957};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task526->add_dep(task527);
  task527->add_dep(task83);
  residualq->add_task(task527);

  auto tensor528 = vector<shared_ptr<Tensor>>{I957, v2_};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task83);
  residualq->add_task(task528);

  vector<IndexRange> I959_index = {virt_, closed_, virt_, active_};
  auto I959 = make_shared<Tensor>(I959_index);
  auto tensor529 = vector<shared_ptr<Tensor>>{I108, t2, I959};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task426->add_dep(task529);
  task529->add_dep(task83);
  residualq->add_task(task529);

  vector<IndexRange> I960_index = {virt_, active_, closed_, virt_};
  auto I960 = make_shared<Tensor>(I960_index);
  auto tensor530 = vector<shared_ptr<Tensor>>{I959, Gamma27_(), I960};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task529->add_dep(task530);
  task530->add_dep(task83);
  residualq->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I960, v2_};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task530->add_dep(task531);
  task531->add_dep(task83);
  residualq->add_task(task531);

  vector<IndexRange> I962_index = {virt_, closed_, virt_, active_};
  auto I962 = make_shared<Tensor>(I962_index);
  auto tensor532 = vector<shared_ptr<Tensor>>{I108, t2, I962};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task426->add_dep(task532);
  task532->add_dep(task83);
  residualq->add_task(task532);

  vector<IndexRange> I963_index = {virt_, active_, closed_, virt_};
  auto I963 = make_shared<Tensor>(I963_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{I962, Gamma27_(), I963};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task532->add_dep(task533);
  task533->add_dep(task83);
  residualq->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I963, v2_};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task83);
  residualq->add_task(task534);
}

#endif
