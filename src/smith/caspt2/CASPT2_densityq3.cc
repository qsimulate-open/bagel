//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityq3.cc
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
#include <src/smith/caspt2/CASPT2_tasks9.h>
#include <src/smith/caspt2/CASPT2_tasks10.h>
#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::make_densityq3(shared_ptr<Queue> densityq, shared_ptr<Task> task285, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I486_index = {active_, virt_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor425 = vector<shared_ptr<Tensor>>{den2, I486};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task425->add_dep(task285);
  densityq->add_task(task425);

  vector<IndexRange> I487_index = {closed_, active_, active_, active_};
  auto I487 = make_shared<Tensor>(I487_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{I486, t2, I487};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  task426->add_dep(task285);
  densityq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I487, Gamma6_(), t2};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task285);
  densityq->add_task(task427);

  vector<IndexRange> I639_index = {virt_, active_, active_, active_};
  auto I639 = make_shared<Tensor>(I639_index);
  auto tensor428 = vector<shared_ptr<Tensor>>{I486, t2, I639};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task425->add_dep(task428);
  task428->add_dep(task285);
  densityq->add_task(task428);

  auto tensor429 = vector<shared_ptr<Tensor>>{I639, Gamma59_(), t2};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task428->add_dep(task429);
  task429->add_dep(task285);
  densityq->add_task(task429);

  vector<IndexRange> I498_index = {closed_, closed_};
  auto I498 = make_shared<Tensor>(I498_index);
  auto tensor430 = vector<shared_ptr<Tensor>>{den2, I498};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task430->add_dep(task285);
  densityq->add_task(task430);

  vector<IndexRange> I499_index = {virt_, closed_, active_, active_};
  auto I499 = make_shared<Tensor>(I499_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{I498, t2, I499};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task430->add_dep(task431);
  task431->add_dep(task285);
  densityq->add_task(task431);

  auto tensor432 = vector<shared_ptr<Tensor>>{I499, Gamma35_(), t2};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task285);
  densityq->add_task(task432);

  vector<IndexRange> I508_index = {virt_, closed_, active_, active_};
  auto I508 = make_shared<Tensor>(I508_index);
  auto tensor433 = vector<shared_ptr<Tensor>>{I498, t2, I508};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task430->add_dep(task433);
  task433->add_dep(task285);
  densityq->add_task(task433);

  auto tensor434 = vector<shared_ptr<Tensor>>{I508, Gamma35_(), t2};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task433->add_dep(task434);
  task434->add_dep(task285);
  densityq->add_task(task434);

  vector<IndexRange> I501_index = {virt_, virt_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{den2, I501};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task435->add_dep(task285);
  densityq->add_task(task435);

  vector<IndexRange> I502_index = {virt_, closed_, active_, active_};
  auto I502 = make_shared<Tensor>(I502_index);
  auto tensor436 = vector<shared_ptr<Tensor>>{I501, t2, I502};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task285);
  densityq->add_task(task436);

  auto tensor437 = vector<shared_ptr<Tensor>>{I502, Gamma35_(), t2};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task436->add_dep(task437);
  task437->add_dep(task285);
  densityq->add_task(task437);

  vector<IndexRange> I511_index = {virt_, closed_, active_, active_};
  auto I511 = make_shared<Tensor>(I511_index);
  auto tensor438 = vector<shared_ptr<Tensor>>{I501, t2, I511};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task435->add_dep(task438);
  task438->add_dep(task285);
  densityq->add_task(task438);

  auto tensor439 = vector<shared_ptr<Tensor>>{I511, Gamma35_(), t2};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task438->add_dep(task439);
  task439->add_dep(task285);
  densityq->add_task(task439);

  vector<IndexRange> I648_index = {virt_, virt_, active_, active_};
  auto I648 = make_shared<Tensor>(I648_index);
  auto tensor440 = vector<shared_ptr<Tensor>>{I501, t2, I648};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task435->add_dep(task440);
  task440->add_dep(task285);
  densityq->add_task(task440);

  auto tensor441 = vector<shared_ptr<Tensor>>{I648, Gamma60_(), t2};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task440->add_dep(task441);
  task441->add_dep(task285);
  densityq->add_task(task441);

  vector<IndexRange> I513_index = {active_, closed_};
  auto I513 = make_shared<Tensor>(I513_index);
  auto tensor442 = vector<shared_ptr<Tensor>>{den2, I513};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task442->add_dep(task285);
  densityq->add_task(task442);

  vector<IndexRange> I514_index = {virt_, active_, active_, active_};
  auto I514 = make_shared<Tensor>(I514_index);
  auto tensor443 = vector<shared_ptr<Tensor>>{I513, t2, I514};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task442->add_dep(task443);
  task443->add_dep(task285);
  densityq->add_task(task443);

  auto tensor444 = vector<shared_ptr<Tensor>>{I514, Gamma51_(), t2};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task443->add_dep(task444);
  task444->add_dep(task285);
  densityq->add_task(task444);

  vector<IndexRange> I537_index = {virt_, virt_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{den2, I537};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task445->add_dep(task285);
  densityq->add_task(task445);

  vector<IndexRange> I538_index = {virt_, active_, active_, active_};
  auto I538 = make_shared<Tensor>(I538_index);
  auto tensor446 = vector<shared_ptr<Tensor>>{I537, t2, I538};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task445->add_dep(task446);
  task446->add_dep(task285);
  densityq->add_task(task446);

  auto tensor447 = vector<shared_ptr<Tensor>>{I538, Gamma59_(), t2};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task446->add_dep(task447);
  task447->add_dep(task285);
  densityq->add_task(task447);

  vector<IndexRange> I549_index = {virt_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor448 = vector<shared_ptr<Tensor>>{den2, I549};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task448->add_dep(task285);
  densityq->add_task(task448);

  vector<IndexRange> I550_index = {virt_, active_};
  auto I550 = make_shared<Tensor>(I550_index);
  auto tensor449 = vector<shared_ptr<Tensor>>{I549, Gamma16_(), I550};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task448->add_dep(task449);
  task449->add_dep(task285);
  densityq->add_task(task449);

  auto tensor450 = vector<shared_ptr<Tensor>>{I550, t2};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task449->add_dep(task450);
  task450->add_dep(task285);
  densityq->add_task(task450);

  vector<IndexRange> I552_index = {virt_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor451 = vector<shared_ptr<Tensor>>{den2, I552};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task451->add_dep(task285);
  densityq->add_task(task451);

  vector<IndexRange> I553_index = {virt_, active_};
  auto I553 = make_shared<Tensor>(I553_index);
  auto tensor452 = vector<shared_ptr<Tensor>>{I552, Gamma16_(), I553};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task451->add_dep(task452);
  task452->add_dep(task285);
  densityq->add_task(task452);

  auto tensor453 = vector<shared_ptr<Tensor>>{I553, t2};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task452->add_dep(task453);
  task453->add_dep(task285);
  densityq->add_task(task453);

  vector<IndexRange> I558_index = {virt_, closed_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{den2, I558};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task454->add_dep(task285);
  densityq->add_task(task454);

  vector<IndexRange> I559_index = {virt_, closed_};
  auto I559 = make_shared<Tensor>(I559_index);
  auto tensor455 = vector<shared_ptr<Tensor>>{I558, t2, I559};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task285);
  densityq->add_task(task455);

  vector<IndexRange> I560_index = {active_, virt_, closed_, active_};
  auto I560 = make_shared<Tensor>(I560_index);
  auto tensor456 = vector<shared_ptr<Tensor>>{I559, Gamma38_(), I560};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task455->add_dep(task456);
  task456->add_dep(task285);
  densityq->add_task(task456);

  auto tensor457 = vector<shared_ptr<Tensor>>{I560, t2};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task456->add_dep(task457);
  task457->add_dep(task285);
  densityq->add_task(task457);

  vector<IndexRange> I567_index = {active_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  auto tensor458 = vector<shared_ptr<Tensor>>{den2, I567};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task458->add_dep(task285);
  densityq->add_task(task458);

  vector<IndexRange> I568_index;
  auto I568 = make_shared<Tensor>(I568_index);
  auto tensor459 = vector<shared_ptr<Tensor>>{I567, Gamma38_(), I568};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task458->add_dep(task459);
  task459->add_dep(task285);
  densityq->add_task(task459);

  auto tensor460 = vector<shared_ptr<Tensor>>{I568, t2};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task459->add_dep(task460);
  task460->add_dep(task285);
  densityq->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I568, t2};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task459->add_dep(task461);
  task461->add_dep(task285);
  densityq->add_task(task461);

  shared_ptr<Tensor> I573;
  if (diagonal) {
    vector<IndexRange> I573_index = {closed_, closed_};
    I573 = make_shared<Tensor>(I573_index);
  }
  shared_ptr<Task462> task462;
  if (diagonal) {
    auto tensor462 = vector<shared_ptr<Tensor>>{den2, I573};
    task462 = make_shared<Task462>(tensor462, pindex);
    task462->add_dep(task285);
    densityq->add_task(task462);
  }

  shared_ptr<Task463> task463;
  if (diagonal) {
    auto tensor463 = vector<shared_ptr<Tensor>>{I573, t2};
    task463 = make_shared<Task463>(tensor463, pindex);
    task462->add_dep(task463);
    task463->add_dep(task285);
    densityq->add_task(task463);
  }

  shared_ptr<Task464> task464;
  if (diagonal) {
    auto tensor464 = vector<shared_ptr<Tensor>>{I573, t2};
    task464 = make_shared<Task464>(tensor464, pindex);
    task462->add_dep(task464);
    task464->add_dep(task285);
    densityq->add_task(task464);
  }

  shared_ptr<Tensor> I577;
  if (diagonal) {
    vector<IndexRange> I577_index = {virt_, virt_};
    I577 = make_shared<Tensor>(I577_index);
  }
  shared_ptr<Task465> task465;
  if (diagonal) {
    auto tensor465 = vector<shared_ptr<Tensor>>{den2, I577};
    task465 = make_shared<Task465>(tensor465, pindex);
    task465->add_dep(task285);
    densityq->add_task(task465);
  }

  shared_ptr<Task466> task466;
  if (diagonal) {
    auto tensor466 = vector<shared_ptr<Tensor>>{I577, t2};
    task466 = make_shared<Task466>(tensor466, pindex);
    task465->add_dep(task466);
    task466->add_dep(task285);
    densityq->add_task(task466);
  }

  shared_ptr<Task467> task467;
  if (diagonal) {
    auto tensor467 = vector<shared_ptr<Tensor>>{I577, t2};
    task467 = make_shared<Task467>(tensor467, pindex);
    task465->add_dep(task467);
    task467->add_dep(task285);
    densityq->add_task(task467);
  }

  vector<IndexRange> I581_index = {closed_, active_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{den2, I581};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task468->add_dep(task285);
  densityq->add_task(task468);

  vector<IndexRange> I582_index = {closed_, active_};
  auto I582 = make_shared<Tensor>(I582_index);
  auto tensor469 = vector<shared_ptr<Tensor>>{I581, Gamma38_(), I582};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task285);
  densityq->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task469->add_dep(task470);
  task470->add_dep(task285);
  densityq->add_task(task470);

  auto tensor471 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task469->add_dep(task471);
  task471->add_dep(task285);
  densityq->add_task(task471);

  vector<IndexRange> I587_index = {active_, virt_};
  auto I587 = make_shared<Tensor>(I587_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{den2, I587};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task472->add_dep(task285);
  densityq->add_task(task472);

  vector<IndexRange> I588_index = {virt_, closed_, active_, active_};
  auto I588 = make_shared<Tensor>(I588_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{I587, t2, I588};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task472->add_dep(task473);
  task473->add_dep(task285);
  densityq->add_task(task473);

  vector<IndexRange> I589_index = {active_, virt_, closed_, active_};
  auto I589 = make_shared<Tensor>(I589_index);
  auto tensor474 = vector<shared_ptr<Tensor>>{I588, Gamma35_(), I589};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task473->add_dep(task474);
  task474->add_dep(task285);
  densityq->add_task(task474);

  auto tensor475 = vector<shared_ptr<Tensor>>{I589, t2};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task474->add_dep(task475);
  task475->add_dep(task285);
  densityq->add_task(task475);

  vector<IndexRange> I590_index = {active_, virt_};
  auto I590 = make_shared<Tensor>(I590_index);
  auto tensor476 = vector<shared_ptr<Tensor>>{den2, I590};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task476->add_dep(task285);
  densityq->add_task(task476);

  vector<IndexRange> I591_index = {virt_, closed_, active_, active_};
  auto I591 = make_shared<Tensor>(I591_index);
  auto tensor477 = vector<shared_ptr<Tensor>>{I590, t2, I591};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task476->add_dep(task477);
  task477->add_dep(task285);
  densityq->add_task(task477);

  auto tensor478 = vector<shared_ptr<Tensor>>{I591, Gamma32_(), t2};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task477->add_dep(task478);
  task478->add_dep(task285);
  densityq->add_task(task478);

  auto tensor479 = vector<shared_ptr<Tensor>>{I591, Gamma35_(), t2};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task477->add_dep(task479);
  task479->add_dep(task285);
  densityq->add_task(task479);

  vector<IndexRange> I599_index = {closed_, virt_};
  auto I599 = make_shared<Tensor>(I599_index);
  auto tensor480 = vector<shared_ptr<Tensor>>{den2, I599};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task480->add_dep(task285);
  densityq->add_task(task480);

  vector<IndexRange> I600_index = {virt_, active_};
  auto I600 = make_shared<Tensor>(I600_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{I599, t2, I600};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task480->add_dep(task481);
  task481->add_dep(task285);
  densityq->add_task(task481);

  auto tensor482 = vector<shared_ptr<Tensor>>{I600, Gamma60_(), t2};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task481->add_dep(task482);
  task482->add_dep(task285);
  densityq->add_task(task482);

  vector<IndexRange> I602_index = {virt_, closed_};
  auto I602 = make_shared<Tensor>(I602_index);
  auto tensor483 = vector<shared_ptr<Tensor>>{den2, I602};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task483->add_dep(task285);
  densityq->add_task(task483);

  vector<IndexRange> I603_index = {virt_, active_};
  auto I603 = make_shared<Tensor>(I603_index);
  auto tensor484 = vector<shared_ptr<Tensor>>{I602, t2, I603};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task483->add_dep(task484);
  task484->add_dep(task285);
  densityq->add_task(task484);

  auto tensor485 = vector<shared_ptr<Tensor>>{I603, Gamma60_(), t2};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task484->add_dep(task485);
  task485->add_dep(task285);
  densityq->add_task(task485);

  vector<IndexRange> I605_index = {closed_, active_};
  auto I605 = make_shared<Tensor>(I605_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{den2, I605};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task486->add_dep(task285);
  densityq->add_task(task486);

  vector<IndexRange> I606_index = {active_, closed_};
  auto I606 = make_shared<Tensor>(I606_index);
  auto tensor487 = vector<shared_ptr<Tensor>>{I605, Gamma38_(), I606};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task285);
  densityq->add_task(task487);

  auto tensor488 = vector<shared_ptr<Tensor>>{I606, t2};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task487->add_dep(task488);
  task488->add_dep(task285);
  densityq->add_task(task488);

  auto tensor489 = vector<shared_ptr<Tensor>>{I606, t2};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task487->add_dep(task489);
  task489->add_dep(task285);
  densityq->add_task(task489);

  vector<IndexRange> I617_index = {closed_, closed_};
  auto I617 = make_shared<Tensor>(I617_index);
  auto tensor490 = vector<shared_ptr<Tensor>>{den2, I617};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task490->add_dep(task285);
  densityq->add_task(task490);

  vector<IndexRange> I618_index = {virt_, closed_, virt_, active_};
  auto I618 = make_shared<Tensor>(I618_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{I617, t2, I618};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task490->add_dep(task491);
  task491->add_dep(task285);
  densityq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I618, Gamma38_(), t2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task285);
  densityq->add_task(task492);

  vector<IndexRange> I621_index = {virt_, closed_, virt_, active_};
  auto I621 = make_shared<Tensor>(I621_index);
  auto tensor493 = vector<shared_ptr<Tensor>>{I617, t2, I621};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task490->add_dep(task493);
  task493->add_dep(task285);
  densityq->add_task(task493);

  auto tensor494 = vector<shared_ptr<Tensor>>{I621, Gamma38_(), t2};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task493->add_dep(task494);
  task494->add_dep(task285);
  densityq->add_task(task494);

  vector<IndexRange> I623_index = {virt_, virt_};
  auto I623 = make_shared<Tensor>(I623_index);
  auto tensor495 = vector<shared_ptr<Tensor>>{den2, I623};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task495->add_dep(task285);
  densityq->add_task(task495);

  vector<IndexRange> I624_index = {virt_, closed_, virt_, active_};
  auto I624 = make_shared<Tensor>(I624_index);
  auto tensor496 = vector<shared_ptr<Tensor>>{I623, t2, I624};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task495->add_dep(task496);
  task496->add_dep(task285);
  densityq->add_task(task496);

  auto tensor497 = vector<shared_ptr<Tensor>>{I624, Gamma38_(), t2};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task496->add_dep(task497);
  task497->add_dep(task285);
  densityq->add_task(task497);

  vector<IndexRange> I627_index = {virt_, closed_, virt_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  auto tensor498 = vector<shared_ptr<Tensor>>{I623, t2, I627};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task495->add_dep(task498);
  task498->add_dep(task285);
  densityq->add_task(task498);

  auto tensor499 = vector<shared_ptr<Tensor>>{I627, Gamma38_(), t2};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task498->add_dep(task499);
  task499->add_dep(task285);
  densityq->add_task(task499);

  vector<IndexRange> I629_index = {virt_, virt_};
  auto I629 = make_shared<Tensor>(I629_index);
  auto tensor500 = vector<shared_ptr<Tensor>>{den2, I629};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task500->add_dep(task285);
  densityq->add_task(task500);

  vector<IndexRange> I630_index = {virt_, closed_, virt_, active_};
  auto I630 = make_shared<Tensor>(I630_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{I629, t2, I630};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task500->add_dep(task501);
  task501->add_dep(task285);
  densityq->add_task(task501);

  auto tensor502 = vector<shared_ptr<Tensor>>{I630, Gamma38_(), t2};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task285);
  densityq->add_task(task502);

  vector<IndexRange> I633_index = {virt_, closed_, virt_, active_};
  auto I633 = make_shared<Tensor>(I633_index);
  auto tensor503 = vector<shared_ptr<Tensor>>{I629, t2, I633};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task500->add_dep(task503);
  task503->add_dep(task285);
  densityq->add_task(task503);

  auto tensor504 = vector<shared_ptr<Tensor>>{I633, Gamma38_(), t2};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task503->add_dep(task504);
  task504->add_dep(task285);
  densityq->add_task(task504);

  vector<IndexRange> I635_index = {active_, closed_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor505 = vector<shared_ptr<Tensor>>{den2, I635};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task505->add_dep(task285);
  densityq->add_task(task505);

  vector<IndexRange> I636_index = {virt_, virt_, active_, active_};
  auto I636 = make_shared<Tensor>(I636_index);
  auto tensor506 = vector<shared_ptr<Tensor>>{I635, t2, I636};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task505->add_dep(task506);
  task506->add_dep(task285);
  densityq->add_task(task506);

  auto tensor507 = vector<shared_ptr<Tensor>>{I636, Gamma60_(), t2};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task506->add_dep(task507);
  task507->add_dep(task285);
  densityq->add_task(task507);
}

#endif
