//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_densityqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_densityq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor428 = {den2};
  auto task428 = make_shared<Task428>(tensor428);
  densityq->add_task(task428);

  vector<IndexRange> I458_index = {active_, active_};
  auto I458 = make_shared<Tensor>(I458_index);
  vector<shared_ptr<Tensor>> tensor429 = {den2, I458};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task429->add_dep(task428);
  densityq->add_task(task429);

  vector<IndexRange> I459_index = {active_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  vector<shared_ptr<Tensor>> tensor430 = {I458, Gamma162_(), I459};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task428);
  densityq->add_task(task430);

  vector<IndexRange> I460_index = {active_, closed_, active_, closed_};
  auto I460 = make_shared<Tensor>(I460_index);
  vector<shared_ptr<Tensor>> tensor431 = {I459, t2, I460};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task430->add_dep(task431);
  task431->add_dep(task428);
  densityq->add_task(task431);

  vector<shared_ptr<Tensor>> tensor432 = {I460, t2};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task428);
  densityq->add_task(task432);

  vector<IndexRange> I552_index = {active_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  vector<shared_ptr<Tensor>> tensor433 = {I458, Gamma193_(), I552};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task429->add_dep(task433);
  task433->add_dep(task428);
  densityq->add_task(task433);

  vector<IndexRange> I553_index = {active_, closed_, virt_, active_};
  auto I553 = make_shared<Tensor>(I553_index);
  vector<shared_ptr<Tensor>> tensor434 = {I552, t2, I553};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task433->add_dep(task434);
  task434->add_dep(task428);
  densityq->add_task(task434);

  vector<shared_ptr<Tensor>> tensor435 = {I553, t2};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task434->add_dep(task435);
  task435->add_dep(task428);
  densityq->add_task(task435);

  vector<IndexRange> I561_index = {active_, active_, active_, active_};
  auto I561 = make_shared<Tensor>(I561_index);
  vector<shared_ptr<Tensor>> tensor436 = {I458, Gamma196_(), I561};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task429->add_dep(task436);
  task436->add_dep(task428);
  densityq->add_task(task436);

  vector<IndexRange> I562_index = {active_, closed_, virt_, active_};
  auto I562 = make_shared<Tensor>(I562_index);
  vector<shared_ptr<Tensor>> tensor437 = {I561, t2, I562};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task436->add_dep(task437);
  task437->add_dep(task428);
  densityq->add_task(task437);

  vector<shared_ptr<Tensor>> tensor438 = {I562, t2};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task428);
  densityq->add_task(task438);

  vector<IndexRange> I595_index = {active_, active_, virt_, closed_};
  auto I595 = make_shared<Tensor>(I595_index);
  vector<shared_ptr<Tensor>> tensor439 = {I561, t2, I595};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task436->add_dep(task439);
  task439->add_dep(task428);
  densityq->add_task(task439);

  vector<shared_ptr<Tensor>> tensor440 = {I595, t2};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task428);
  densityq->add_task(task440);

  vector<IndexRange> I604_index = {active_, active_, virt_, closed_};
  auto I604 = make_shared<Tensor>(I604_index);
  vector<shared_ptr<Tensor>> tensor441 = {I561, t2, I604};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task436->add_dep(task441);
  task441->add_dep(task428);
  densityq->add_task(task441);

  vector<shared_ptr<Tensor>> tensor442 = {I604, t2};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task441->add_dep(task442);
  task442->add_dep(task428);
  densityq->add_task(task442);

  vector<IndexRange> I743_index = {active_, active_, active_, active_};
  auto I743 = make_shared<Tensor>(I743_index);
  vector<shared_ptr<Tensor>> tensor443 = {I458, Gamma254_(), I743};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task429->add_dep(task443);
  task443->add_dep(task428);
  densityq->add_task(task443);

  vector<IndexRange> I744_index = {virt_, active_, virt_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  vector<shared_ptr<Tensor>> tensor444 = {I743, t2, I744};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task443->add_dep(task444);
  task444->add_dep(task428);
  densityq->add_task(task444);

  vector<shared_ptr<Tensor>> tensor445 = {I744, t2};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task444->add_dep(task445);
  task445->add_dep(task428);
  densityq->add_task(task445);

  vector<IndexRange> I461_index = {closed_, closed_};
  auto I461 = make_shared<Tensor>(I461_index);
  vector<shared_ptr<Tensor>> tensor446 = {den2, I461};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task446->add_dep(task428);
  densityq->add_task(task446);

  vector<IndexRange> I462_index = {closed_, closed_, active_, active_};
  auto I462 = make_shared<Tensor>(I462_index);
  vector<shared_ptr<Tensor>> tensor447 = {I461, t2, I462};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task446->add_dep(task447);
  task447->add_dep(task428);
  densityq->add_task(task447);

  vector<IndexRange> I463_index = {active_, closed_, active_, closed_};
  auto I463 = make_shared<Tensor>(I463_index);
  vector<shared_ptr<Tensor>> tensor448 = {I462, Gamma94_(), I463};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task428);
  densityq->add_task(task448);

  vector<shared_ptr<Tensor>> tensor449 = {I463, t2};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task448->add_dep(task449);
  task449->add_dep(task428);
  densityq->add_task(task449);

  vector<IndexRange> I555_index = {closed_, virt_, active_, active_};
  auto I555 = make_shared<Tensor>(I555_index);
  vector<shared_ptr<Tensor>> tensor450 = {I461, t2, I555};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task446->add_dep(task450);
  task450->add_dep(task428);
  densityq->add_task(task450);

  vector<IndexRange> I556_index = {active_, closed_, virt_, active_};
  auto I556 = make_shared<Tensor>(I556_index);
  vector<shared_ptr<Tensor>> tensor451 = {I555, Gamma32_(), I556};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task428);
  densityq->add_task(task451);

  vector<shared_ptr<Tensor>> tensor452 = {I556, t2};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task451->add_dep(task452);
  task452->add_dep(task428);
  densityq->add_task(task452);

  vector<IndexRange> I564_index = {closed_, virt_, active_, active_};
  auto I564 = make_shared<Tensor>(I564_index);
  vector<shared_ptr<Tensor>> tensor453 = {I461, t2, I564};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task446->add_dep(task453);
  task453->add_dep(task428);
  densityq->add_task(task453);

  vector<IndexRange> I565_index = {active_, closed_, virt_, active_};
  auto I565 = make_shared<Tensor>(I565_index);
  vector<shared_ptr<Tensor>> tensor454 = {I564, Gamma35_(), I565};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task453->add_dep(task454);
  task454->add_dep(task428);
  densityq->add_task(task454);

  vector<shared_ptr<Tensor>> tensor455 = {I565, t2};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task428);
  densityq->add_task(task455);

  vector<IndexRange> I464_index = {active_, closed_};
  auto I464 = make_shared<Tensor>(I464_index);
  vector<shared_ptr<Tensor>> tensor456 = {den2, I464};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task456->add_dep(task428);
  densityq->add_task(task456);

  vector<IndexRange> I465_index = {closed_, active_, active_, active_};
  auto I465 = make_shared<Tensor>(I465_index);
  vector<shared_ptr<Tensor>> tensor457 = {I464, t2, I465};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task456->add_dep(task457);
  task457->add_dep(task428);
  densityq->add_task(task457);

  vector<IndexRange> I466_index = {active_, active_, closed_, active_};
  auto I466 = make_shared<Tensor>(I466_index);
  vector<shared_ptr<Tensor>> tensor458 = {I465, Gamma2_(), I466};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task428);
  densityq->add_task(task458);

  vector<shared_ptr<Tensor>> tensor459 = {I466, t2};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task458->add_dep(task459);
  task459->add_dep(task428);
  densityq->add_task(task459);

  vector<IndexRange> I570_index = {virt_, active_, active_, active_};
  auto I570 = make_shared<Tensor>(I570_index);
  vector<shared_ptr<Tensor>> tensor460 = {I464, t2, I570};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task456->add_dep(task460);
  task460->add_dep(task428);
  densityq->add_task(task460);

  vector<IndexRange> I571_index = {active_, virt_, active_, active_};
  auto I571 = make_shared<Tensor>(I571_index);
  vector<shared_ptr<Tensor>> tensor461 = {I570, Gamma37_(), I571};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task460->add_dep(task461);
  task461->add_dep(task428);
  densityq->add_task(task461);

  vector<shared_ptr<Tensor>> tensor462 = {I571, t2};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task461->add_dep(task462);
  task462->add_dep(task428);
  densityq->add_task(task462);

  vector<IndexRange> I467_index = {active_, virt_};
  auto I467 = make_shared<Tensor>(I467_index);
  vector<shared_ptr<Tensor>> tensor463 = {den2, I467};
  auto task463 = make_shared<Task463>(tensor463, pindex);
  task463->add_dep(task428);
  densityq->add_task(task463);

  vector<IndexRange> I468_index = {closed_, closed_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  vector<shared_ptr<Tensor>> tensor464 = {I467, t2, I468};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task463->add_dep(task464);
  task464->add_dep(task428);
  densityq->add_task(task464);

  vector<IndexRange> I469_index = {active_, closed_, active_, closed_};
  auto I469 = make_shared<Tensor>(I469_index);
  vector<shared_ptr<Tensor>> tensor465 = {I468, Gamma3_(), I469};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task464->add_dep(task465);
  task465->add_dep(task428);
  densityq->add_task(task465);

  vector<shared_ptr<Tensor>> tensor466 = {I469, t2};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task465->add_dep(task466);
  task466->add_dep(task428);
  densityq->add_task(task466);

  vector<IndexRange> I579_index = {closed_, virt_, active_, active_};
  auto I579 = make_shared<Tensor>(I579_index);
  vector<shared_ptr<Tensor>> tensor467 = {I467, t2, I579};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task463->add_dep(task467);
  task467->add_dep(task428);
  densityq->add_task(task467);

  vector<IndexRange> I580_index = {active_, closed_, virt_, active_};
  auto I580 = make_shared<Tensor>(I580_index);
  vector<shared_ptr<Tensor>> tensor468 = {I579, Gamma35_(), I580};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task467->add_dep(task468);
  task468->add_dep(task428);
  densityq->add_task(task468);

  vector<shared_ptr<Tensor>> tensor469 = {I580, t2};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task428);
  densityq->add_task(task469);

  vector<IndexRange> I582_index = {closed_, virt_, active_, active_};
  auto I582 = make_shared<Tensor>(I582_index);
  vector<shared_ptr<Tensor>> tensor470 = {I467, t2, I582};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task463->add_dep(task470);
  task470->add_dep(task428);
  densityq->add_task(task470);

  vector<IndexRange> I583_index = {active_, closed_, virt_, active_};
  auto I583 = make_shared<Tensor>(I583_index);
  vector<shared_ptr<Tensor>> tensor471 = {I582, Gamma32_(), I583};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task470->add_dep(task471);
  task471->add_dep(task428);
  densityq->add_task(task471);

  vector<shared_ptr<Tensor>> tensor472 = {I583, t2};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task471->add_dep(task472);
  task472->add_dep(task428);
  densityq->add_task(task472);

  vector<IndexRange> I621_index = {virt_, closed_, active_, active_};
  auto I621 = make_shared<Tensor>(I621_index);
  vector<shared_ptr<Tensor>> tensor473 = {I467, t2, I621};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task463->add_dep(task473);
  task473->add_dep(task428);
  densityq->add_task(task473);

  vector<IndexRange> I622_index = {active_, active_, virt_, closed_};
  auto I622 = make_shared<Tensor>(I622_index);
  vector<shared_ptr<Tensor>> tensor474 = {I621, Gamma35_(), I622};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task473->add_dep(task474);
  task474->add_dep(task428);
  densityq->add_task(task474);

  vector<shared_ptr<Tensor>> tensor475 = {I622, t2};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task474->add_dep(task475);
  task475->add_dep(task428);
  densityq->add_task(task475);

  vector<IndexRange> I624_index = {virt_, closed_, active_, active_};
  auto I624 = make_shared<Tensor>(I624_index);
  vector<shared_ptr<Tensor>> tensor476 = {I467, t2, I624};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task463->add_dep(task476);
  task476->add_dep(task428);
  densityq->add_task(task476);

  vector<IndexRange> I625_index = {active_, active_, virt_, closed_};
  auto I625 = make_shared<Tensor>(I625_index);
  vector<shared_ptr<Tensor>> tensor477 = {I624, Gamma35_(), I625};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task476->add_dep(task477);
  task477->add_dep(task428);
  densityq->add_task(task477);

  vector<shared_ptr<Tensor>> tensor478 = {I625, t2};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task477->add_dep(task478);
  task478->add_dep(task428);
  densityq->add_task(task478);

  vector<IndexRange> I470_index = {active_, closed_};
  auto I470 = make_shared<Tensor>(I470_index);
  vector<shared_ptr<Tensor>> tensor479 = {den2, I470};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task479->add_dep(task428);
  densityq->add_task(task479);

  vector<IndexRange> I471_index = {closed_, active_, active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  vector<shared_ptr<Tensor>> tensor480 = {I470, t2, I471};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task479->add_dep(task480);
  task480->add_dep(task428);
  densityq->add_task(task480);

  vector<IndexRange> I472_index = {active_, closed_, active_, active_};
  auto I472 = make_shared<Tensor>(I472_index);
  vector<shared_ptr<Tensor>> tensor481 = {I471, Gamma4_(), I472};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task480->add_dep(task481);
  task481->add_dep(task428);
  densityq->add_task(task481);

  vector<shared_ptr<Tensor>> tensor482 = {I472, t2};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task481->add_dep(task482);
  task482->add_dep(task428);
  densityq->add_task(task482);

  vector<IndexRange> I627_index = {virt_, active_, active_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  vector<shared_ptr<Tensor>> tensor483 = {I470, t2, I627};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task479->add_dep(task483);
  task483->add_dep(task428);
  densityq->add_task(task483);

  vector<IndexRange> I628_index = {active_, active_, virt_, active_};
  auto I628 = make_shared<Tensor>(I628_index);
  vector<shared_ptr<Tensor>> tensor484 = {I627, Gamma56_(), I628};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task483->add_dep(task484);
  task484->add_dep(task428);
  densityq->add_task(task484);

  vector<shared_ptr<Tensor>> tensor485 = {I628, t2};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task484->add_dep(task485);
  task485->add_dep(task428);
  densityq->add_task(task485);

  vector<IndexRange> I630_index = {virt_, active_, active_, active_};
  auto I630 = make_shared<Tensor>(I630_index);
  vector<shared_ptr<Tensor>> tensor486 = {I470, t2, I630};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task479->add_dep(task486);
  task486->add_dep(task428);
  densityq->add_task(task486);

  vector<IndexRange> I631_index = {active_, active_, virt_, active_};
  auto I631 = make_shared<Tensor>(I631_index);
  vector<shared_ptr<Tensor>> tensor487 = {I630, Gamma57_(), I631};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task428);
  densityq->add_task(task487);

  vector<shared_ptr<Tensor>> tensor488 = {I631, t2};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task487->add_dep(task488);
  task488->add_dep(task428);
  densityq->add_task(task488);

  vector<IndexRange> I473_index = {active_, active_};
  auto I473 = make_shared<Tensor>(I473_index);
  vector<shared_ptr<Tensor>> tensor489 = {den2, I473};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task489->add_dep(task428);
  densityq->add_task(task489);

  vector<IndexRange> I474_index = {active_, active_, active_, active_, active_, active_};
  auto I474 = make_shared<Tensor>(I474_index);
  vector<shared_ptr<Tensor>> tensor490 = {I473, Gamma167_(), I474};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task489->add_dep(task490);
  task490->add_dep(task428);
  densityq->add_task(task490);

  vector<IndexRange> I475_index = {active_, closed_, active_, active_};
  auto I475 = make_shared<Tensor>(I475_index);
  vector<shared_ptr<Tensor>> tensor491 = {I474, t2, I475};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task490->add_dep(task491);
  task491->add_dep(task428);
  densityq->add_task(task491);

  vector<shared_ptr<Tensor>> tensor492 = {I475, t2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task428);
  densityq->add_task(task492);

  vector<IndexRange> I633_index = {active_, active_, active_, active_, active_, active_};
  auto I633 = make_shared<Tensor>(I633_index);
  vector<shared_ptr<Tensor>> tensor493 = {I473, Gamma220_(), I633};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task489->add_dep(task493);
  task493->add_dep(task428);
  densityq->add_task(task493);

  vector<IndexRange> I634_index = {active_, active_, virt_, active_};
  auto I634 = make_shared<Tensor>(I634_index);
  vector<shared_ptr<Tensor>> tensor494 = {I633, t2, I634};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task493->add_dep(task494);
  task494->add_dep(task428);
  densityq->add_task(task494);

  vector<shared_ptr<Tensor>> tensor495 = {I634, t2};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task494->add_dep(task495);
  task495->add_dep(task428);
  densityq->add_task(task495);

  vector<IndexRange> I476_index = {closed_, closed_};
  auto I476 = make_shared<Tensor>(I476_index);
  vector<shared_ptr<Tensor>> tensor496 = {den2, I476};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task496->add_dep(task428);
  densityq->add_task(task496);

  vector<IndexRange> I477_index = {closed_, active_, active_, active_};
  auto I477 = make_shared<Tensor>(I477_index);
  vector<shared_ptr<Tensor>> tensor497 = {I476, t2, I477};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task496->add_dep(task497);
  task497->add_dep(task428);
  densityq->add_task(task497);

  vector<IndexRange> I478_index = {active_, closed_, active_, active_};
  auto I478 = make_shared<Tensor>(I478_index);
  vector<shared_ptr<Tensor>> tensor498 = {I477, Gamma6_(), I478};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task428);
  densityq->add_task(task498);

  vector<shared_ptr<Tensor>> tensor499 = {I478, t2};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task498->add_dep(task499);
  task499->add_dep(task428);
  densityq->add_task(task499);

  vector<IndexRange> I479_index = {closed_, virt_};
  auto I479 = make_shared<Tensor>(I479_index);
  vector<shared_ptr<Tensor>> tensor500 = {den2, I479};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task500->add_dep(task428);
  densityq->add_task(task500);

  vector<IndexRange> I480_index = {closed_, active_};
  auto I480 = make_shared<Tensor>(I480_index);
  vector<shared_ptr<Tensor>> tensor501 = {I479, t2, I480};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task500->add_dep(task501);
  task501->add_dep(task428);
  densityq->add_task(task501);

  vector<IndexRange> I481_index = {active_, closed_, active_, active_};
  auto I481 = make_shared<Tensor>(I481_index);
  vector<shared_ptr<Tensor>> tensor502 = {I480, Gamma7_(), I481};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task428);
  densityq->add_task(task502);

  vector<shared_ptr<Tensor>> tensor503 = {I481, t2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task502->add_dep(task503);
  task503->add_dep(task428);
  densityq->add_task(task503);

  vector<IndexRange> I483_index = {closed_, active_};
  auto I483 = make_shared<Tensor>(I483_index);
  vector<shared_ptr<Tensor>> tensor504 = {I479, t2, I483};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task500->add_dep(task504);
  task504->add_dep(task428);
  densityq->add_task(task504);

  vector<IndexRange> I484_index = {active_, closed_, active_, active_};
  auto I484 = make_shared<Tensor>(I484_index);
  vector<shared_ptr<Tensor>> tensor505 = {I483, Gamma7_(), I484};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task428);
  densityq->add_task(task505);

  vector<shared_ptr<Tensor>> tensor506 = {I484, t2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task505->add_dep(task506);
  task506->add_dep(task428);
  densityq->add_task(task506);

  vector<IndexRange> I639_index = {virt_, active_};
  auto I639 = make_shared<Tensor>(I639_index);
  vector<shared_ptr<Tensor>> tensor507 = {I479, t2, I639};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task500->add_dep(task507);
  task507->add_dep(task428);
  densityq->add_task(task507);

  vector<IndexRange> I640_index = {active_, active_, virt_, active_};
  auto I640 = make_shared<Tensor>(I640_index);
  vector<shared_ptr<Tensor>> tensor508 = {I639, Gamma60_(), I640};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task507->add_dep(task508);
  task508->add_dep(task428);
  densityq->add_task(task508);

  vector<shared_ptr<Tensor>> tensor509 = {I640, t2};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task508->add_dep(task509);
  task509->add_dep(task428);
  densityq->add_task(task509);

  vector<IndexRange> I642_index = {virt_, active_};
  auto I642 = make_shared<Tensor>(I642_index);
  vector<shared_ptr<Tensor>> tensor510 = {I479, t2, I642};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task500->add_dep(task510);
  task510->add_dep(task428);
  densityq->add_task(task510);

  vector<IndexRange> I643_index = {active_, active_, virt_, active_};
  auto I643 = make_shared<Tensor>(I643_index);
  vector<shared_ptr<Tensor>> tensor511 = {I642, Gamma60_(), I643};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task510->add_dep(task511);
  task511->add_dep(task428);
  densityq->add_task(task511);

  vector<shared_ptr<Tensor>> tensor512 = {I643, t2};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task511->add_dep(task512);
  task512->add_dep(task428);
  densityq->add_task(task512);

  vector<IndexRange> I485_index = {active_, virt_};
  auto I485 = make_shared<Tensor>(I485_index);
  vector<shared_ptr<Tensor>> tensor513 = {den2, I485};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task513->add_dep(task428);
  densityq->add_task(task513);

  vector<IndexRange> I486_index = {closed_, active_, active_, active_};
  auto I486 = make_shared<Tensor>(I486_index);
  vector<shared_ptr<Tensor>> tensor514 = {I485, t2, I486};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task428);
  densityq->add_task(task514);

  vector<IndexRange> I487_index = {active_, closed_, active_, active_};
  auto I487 = make_shared<Tensor>(I487_index);
  vector<shared_ptr<Tensor>> tensor515 = {I486, Gamma9_(), I487};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task514->add_dep(task515);
  task515->add_dep(task428);
  densityq->add_task(task515);

  vector<shared_ptr<Tensor>> tensor516 = {I487, t2};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task515->add_dep(task516);
  task516->add_dep(task428);
  densityq->add_task(task516);

  vector<IndexRange> I489_index = {closed_, active_, active_, active_};
  auto I489 = make_shared<Tensor>(I489_index);
  vector<shared_ptr<Tensor>> tensor517 = {I485, t2, I489};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task513->add_dep(task517);
  task517->add_dep(task428);
  densityq->add_task(task517);

  vector<IndexRange> I490_index = {active_, closed_, active_, active_};
  auto I490 = make_shared<Tensor>(I490_index);
  vector<shared_ptr<Tensor>> tensor518 = {I489, Gamma6_(), I490};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task517->add_dep(task518);
  task518->add_dep(task428);
  densityq->add_task(task518);

  vector<shared_ptr<Tensor>> tensor519 = {I490, t2};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task518->add_dep(task519);
  task519->add_dep(task428);
  densityq->add_task(task519);

  vector<IndexRange> I645_index = {virt_, active_, active_, active_};
  auto I645 = make_shared<Tensor>(I645_index);
  vector<shared_ptr<Tensor>> tensor520 = {I485, t2, I645};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task513->add_dep(task520);
  task520->add_dep(task428);
  densityq->add_task(task520);

  vector<IndexRange> I646_index = {active_, active_, virt_, active_};
  auto I646 = make_shared<Tensor>(I646_index);
  vector<shared_ptr<Tensor>> tensor521 = {I645, Gamma59_(), I646};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task520->add_dep(task521);
  task521->add_dep(task428);
  densityq->add_task(task521);

  vector<shared_ptr<Tensor>> tensor522 = {I646, t2};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task428);
  densityq->add_task(task522);

  vector<IndexRange> I491_index = {active_, virt_};
  auto I491 = make_shared<Tensor>(I491_index);
  vector<shared_ptr<Tensor>> tensor523 = {den2, I491};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task523->add_dep(task428);
  densityq->add_task(task523);

  vector<IndexRange> I492_index = {closed_, closed_, active_, active_};
  auto I492 = make_shared<Tensor>(I492_index);
  vector<shared_ptr<Tensor>> tensor524 = {I491, t2, I492};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task428);
  densityq->add_task(task524);

  vector<IndexRange> I493_index = {closed_, active_, closed_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  vector<shared_ptr<Tensor>> tensor525 = {I492, Gamma3_(), I493};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task428);
  densityq->add_task(task525);

  vector<shared_ptr<Tensor>> tensor526 = {I493, t2};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task525->add_dep(task526);
  task526->add_dep(task428);
  densityq->add_task(task526);

  vector<IndexRange> I494_index = {virt_, closed_};
  auto I494 = make_shared<Tensor>(I494_index);
  vector<shared_ptr<Tensor>> tensor527 = {den2, I494};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task527->add_dep(task428);
  densityq->add_task(task527);

  vector<IndexRange> I495_index = {closed_, active_};
  auto I495 = make_shared<Tensor>(I495_index);
  vector<shared_ptr<Tensor>> tensor528 = {I494, t2, I495};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task428);
  densityq->add_task(task528);

  vector<IndexRange> I496_index = {active_, active_, closed_, active_};
  auto I496 = make_shared<Tensor>(I496_index);
  vector<shared_ptr<Tensor>> tensor529 = {I495, Gamma12_(), I496};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task528->add_dep(task529);
  task529->add_dep(task428);
  densityq->add_task(task529);

  vector<shared_ptr<Tensor>> tensor530 = {I496, t2};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task529->add_dep(task530);
  task530->add_dep(task428);
  densityq->add_task(task530);

  vector<IndexRange> I497_index = {closed_, virt_};
  auto I497 = make_shared<Tensor>(I497_index);
  vector<shared_ptr<Tensor>> tensor531 = {den2, I497};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task531->add_dep(task428);
  densityq->add_task(task531);

  vector<IndexRange> I498_index = {closed_, active_};
  auto I498 = make_shared<Tensor>(I498_index);
  vector<shared_ptr<Tensor>> tensor532 = {I497, t2, I498};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task531->add_dep(task532);
  task532->add_dep(task428);
  densityq->add_task(task532);

  vector<IndexRange> I499_index = {active_, active_, closed_, active_};
  auto I499 = make_shared<Tensor>(I499_index);
  vector<shared_ptr<Tensor>> tensor533 = {I498, Gamma12_(), I499};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task532->add_dep(task533);
  task533->add_dep(task428);
  densityq->add_task(task533);

  vector<shared_ptr<Tensor>> tensor534 = {I499, t2};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task428);
  densityq->add_task(task534);

  vector<IndexRange> I654_index = {virt_, closed_};
  auto I654 = make_shared<Tensor>(I654_index);
  vector<shared_ptr<Tensor>> tensor535 = {I497, t2, I654};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task531->add_dep(task535);
  task535->add_dep(task428);
  densityq->add_task(task535);

  vector<IndexRange> I655_index = {active_, virt_, closed_, active_};
  auto I655 = make_shared<Tensor>(I655_index);
  vector<shared_ptr<Tensor>> tensor536 = {I654, Gamma38_(), I655};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task428);
  densityq->add_task(task536);

  vector<shared_ptr<Tensor>> tensor537 = {I655, t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task428);
  densityq->add_task(task537);

  vector<IndexRange> I500_index = {active_, active_};
  auto I500 = make_shared<Tensor>(I500_index);
  vector<shared_ptr<Tensor>> tensor538 = {den2, I500};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task538->add_dep(task428);
  densityq->add_task(task538);

  vector<IndexRange> I501_index = {active_, active_};
  auto I501 = make_shared<Tensor>(I501_index);
  vector<shared_ptr<Tensor>> tensor539 = {I500, Gamma176_(), I501};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task428);
  densityq->add_task(task539);

  vector<IndexRange> I502_index = {active_, closed_, virt_, closed_};
  auto I502 = make_shared<Tensor>(I502_index);
  vector<shared_ptr<Tensor>> tensor540 = {I501, t2, I502};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task539->add_dep(task540);
  task540->add_dep(task428);
  densityq->add_task(task540);

  vector<shared_ptr<Tensor>> tensor541 = {I502, t2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task428);
  densityq->add_task(task541);

  vector<IndexRange> I505_index = {active_, closed_, virt_, closed_};
  auto I505 = make_shared<Tensor>(I505_index);
  vector<shared_ptr<Tensor>> tensor542 = {I501, t2, I505};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task539->add_dep(task542);
  task542->add_dep(task428);
  densityq->add_task(task542);

  vector<shared_ptr<Tensor>> tensor543 = {I505, t2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task428);
  densityq->add_task(task543);

  vector<IndexRange> I710_index = {active_, active_};
  auto I710 = make_shared<Tensor>(I710_index);
  vector<shared_ptr<Tensor>> tensor544 = {I500, Gamma60_(), I710};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task538->add_dep(task544);
  task544->add_dep(task428);
  densityq->add_task(task544);

  vector<IndexRange> I711_index = {virt_, closed_, virt_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  vector<shared_ptr<Tensor>> tensor545 = {I710, t2, I711};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task544->add_dep(task545);
  task545->add_dep(task428);
  densityq->add_task(task545);

  vector<shared_ptr<Tensor>> tensor546 = {I711, t2};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task545->add_dep(task546);
  task546->add_dep(task428);
  densityq->add_task(task546);

  vector<IndexRange> I714_index = {virt_, closed_, virt_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  vector<shared_ptr<Tensor>> tensor547 = {I710, t2, I714};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task544->add_dep(task547);
  task547->add_dep(task428);
  densityq->add_task(task547);

  vector<shared_ptr<Tensor>> tensor548 = {I714, t2};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task547->add_dep(task548);
  task548->add_dep(task428);
  densityq->add_task(task548);

  vector<IndexRange> I506_index = {closed_, closed_};
  auto I506 = make_shared<Tensor>(I506_index);
  vector<shared_ptr<Tensor>> tensor549 = {den2, I506};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task549->add_dep(task428);
  densityq->add_task(task549);

  vector<IndexRange> I507_index = {closed_, virt_, closed_, active_};
  auto I507 = make_shared<Tensor>(I507_index);
  vector<shared_ptr<Tensor>> tensor550 = {I506, t2, I507};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task549->add_dep(task550);
  task550->add_dep(task428);
  densityq->add_task(task550);

  vector<IndexRange> I508_index = {active_, closed_, virt_, closed_};
  auto I508 = make_shared<Tensor>(I508_index);
  vector<shared_ptr<Tensor>> tensor551 = {I507, Gamma16_(), I508};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task550->add_dep(task551);
  task551->add_dep(task428);
  densityq->add_task(task551);

  vector<shared_ptr<Tensor>> tensor552 = {I508, t2};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task551->add_dep(task552);
  task552->add_dep(task428);
  densityq->add_task(task552);

  vector<IndexRange> I510_index = {closed_, virt_, closed_, active_};
  auto I510 = make_shared<Tensor>(I510_index);
  vector<shared_ptr<Tensor>> tensor553 = {I506, t2, I510};
  auto task553 = make_shared<Task553>(tensor553, pindex);
  task549->add_dep(task553);
  task553->add_dep(task428);
  densityq->add_task(task553);

  vector<IndexRange> I511_index = {active_, closed_, virt_, closed_};
  auto I511 = make_shared<Tensor>(I511_index);
  vector<shared_ptr<Tensor>> tensor554 = {I510, Gamma16_(), I511};
  auto task554 = make_shared<Task554>(tensor554, pindex);
  task553->add_dep(task554);
  task554->add_dep(task428);
  densityq->add_task(task554);

  vector<shared_ptr<Tensor>> tensor555 = {I511, t2};
  auto task555 = make_shared<Task555>(tensor555, pindex);
  task554->add_dep(task555);
  task555->add_dep(task428);
  densityq->add_task(task555);

  vector<IndexRange> I512_index = {closed_, closed_};
  auto I512 = make_shared<Tensor>(I512_index);
  vector<shared_ptr<Tensor>> tensor556 = {den2, I512};
  auto task556 = make_shared<Task556>(tensor556, pindex);
  task556->add_dep(task428);
  densityq->add_task(task556);

  vector<IndexRange> I513_index = {closed_, virt_, closed_, active_};
  auto I513 = make_shared<Tensor>(I513_index);
  vector<shared_ptr<Tensor>> tensor557 = {I512, t2, I513};
  auto task557 = make_shared<Task557>(tensor557, pindex);
  task556->add_dep(task557);
  task557->add_dep(task428);
  densityq->add_task(task557);

  vector<IndexRange> I514_index = {active_, closed_, virt_, closed_};
  auto I514 = make_shared<Tensor>(I514_index);
  vector<shared_ptr<Tensor>> tensor558 = {I513, Gamma16_(), I514};
  auto task558 = make_shared<Task558>(tensor558, pindex);
  task557->add_dep(task558);
  task558->add_dep(task428);
  densityq->add_task(task558);

  vector<shared_ptr<Tensor>> tensor559 = {I514, t2};
  auto task559 = make_shared<Task559>(tensor559, pindex);
  task558->add_dep(task559);
  task559->add_dep(task428);
  densityq->add_task(task559);

  vector<IndexRange> I519_index = {closed_, virt_, closed_, active_};
  auto I519 = make_shared<Tensor>(I519_index);
  vector<shared_ptr<Tensor>> tensor560 = {I512, t2, I519};
  auto task560 = make_shared<Task560>(tensor560, pindex);
  task556->add_dep(task560);
  task560->add_dep(task428);
  densityq->add_task(task560);

  vector<IndexRange> I520_index = {active_, closed_, virt_, closed_};
  auto I520 = make_shared<Tensor>(I520_index);
  vector<shared_ptr<Tensor>> tensor561 = {I519, Gamma16_(), I520};
  auto task561 = make_shared<Task561>(tensor561, pindex);
  task560->add_dep(task561);
  task561->add_dep(task428);
  densityq->add_task(task561);

  vector<shared_ptr<Tensor>> tensor562 = {I520, t2};
  auto task562 = make_shared<Task562>(tensor562, pindex);
  task561->add_dep(task562);
  task562->add_dep(task428);
  densityq->add_task(task562);

  vector<IndexRange> I515_index = {virt_, virt_};
  auto I515 = make_shared<Tensor>(I515_index);
  vector<shared_ptr<Tensor>> tensor563 = {den2, I515};
  auto task563 = make_shared<Task563>(tensor563, pindex);
  task563->add_dep(task428);
  densityq->add_task(task563);

  vector<IndexRange> I516_index = {closed_, virt_, closed_, active_};
  auto I516 = make_shared<Tensor>(I516_index);
  vector<shared_ptr<Tensor>> tensor564 = {I515, t2, I516};
  auto task564 = make_shared<Task564>(tensor564, pindex);
  task563->add_dep(task564);
  task564->add_dep(task428);
  densityq->add_task(task564);

  vector<IndexRange> I517_index = {active_, closed_, virt_, closed_};
  auto I517 = make_shared<Tensor>(I517_index);
  vector<shared_ptr<Tensor>> tensor565 = {I516, Gamma16_(), I517};
  auto task565 = make_shared<Task565>(tensor565, pindex);
  task564->add_dep(task565);
  task565->add_dep(task428);
  densityq->add_task(task565);

  vector<shared_ptr<Tensor>> tensor566 = {I517, t2};
  auto task566 = make_shared<Task566>(tensor566, pindex);
  task565->add_dep(task566);
  task566->add_dep(task428);
  densityq->add_task(task566);

  vector<IndexRange> I522_index = {closed_, virt_, closed_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  vector<shared_ptr<Tensor>> tensor567 = {I515, t2, I522};
  auto task567 = make_shared<Task567>(tensor567, pindex);
  task563->add_dep(task567);
  task567->add_dep(task428);
  densityq->add_task(task567);

  vector<IndexRange> I523_index = {active_, closed_, virt_, closed_};
  auto I523 = make_shared<Tensor>(I523_index);
  vector<shared_ptr<Tensor>> tensor568 = {I522, Gamma16_(), I523};
  auto task568 = make_shared<Task568>(tensor568, pindex);
  task567->add_dep(task568);
  task568->add_dep(task428);
  densityq->add_task(task568);

  vector<shared_ptr<Tensor>> tensor569 = {I523, t2};
  auto task569 = make_shared<Task569>(tensor569, pindex);
  task568->add_dep(task569);
  task569->add_dep(task428);
  densityq->add_task(task569);

  vector<IndexRange> I524_index = {active_, closed_};
  auto I524 = make_shared<Tensor>(I524_index);
  vector<shared_ptr<Tensor>> tensor570 = {den2, I524};
  auto task570 = make_shared<Task570>(tensor570, pindex);
  task570->add_dep(task428);
  densityq->add_task(task570);

  vector<IndexRange> I525_index = {virt_, closed_, active_, active_};
  auto I525 = make_shared<Tensor>(I525_index);
  vector<shared_ptr<Tensor>> tensor571 = {I524, t2, I525};
  auto task571 = make_shared<Task571>(tensor571, pindex);
  task570->add_dep(task571);
  task571->add_dep(task428);
  densityq->add_task(task571);

  vector<IndexRange> I526_index = {active_, virt_, closed_, active_};
  auto I526 = make_shared<Tensor>(I526_index);
  vector<shared_ptr<Tensor>> tensor572 = {I525, Gamma22_(), I526};
  auto task572 = make_shared<Task572>(tensor572, pindex);
  task571->add_dep(task572);
  task572->add_dep(task428);
  densityq->add_task(task572);

  vector<shared_ptr<Tensor>> tensor573 = {I526, t2};
  auto task573 = make_shared<Task573>(tensor573, pindex);
  task572->add_dep(task573);
  task573->add_dep(task428);
  densityq->add_task(task573);

  vector<IndexRange> I532_index = {closed_, virt_, active_, active_};
  auto I532 = make_shared<Tensor>(I532_index);
  vector<shared_ptr<Tensor>> tensor574 = {I525, Gamma12_(), I532};
  auto task574 = make_shared<Task574>(tensor574, pindex);
  task571->add_dep(task574);
  task574->add_dep(task428);
  densityq->add_task(task574);

  vector<shared_ptr<Tensor>> tensor575 = {I532, t2};
  auto task575 = make_shared<Task575>(tensor575, pindex);
  task574->add_dep(task575);
  task575->add_dep(task428);
  densityq->add_task(task575);

  vector<IndexRange> I527_index = {active_, closed_};
  auto I527 = make_shared<Tensor>(I527_index);
  vector<shared_ptr<Tensor>> tensor576 = {den2, I527};
  auto task576 = make_shared<Task576>(tensor576, pindex);
  task576->add_dep(task428);
  densityq->add_task(task576);

  vector<IndexRange> I528_index = {virt_, closed_, active_, active_};
  auto I528 = make_shared<Tensor>(I528_index);
  vector<shared_ptr<Tensor>> tensor577 = {I527, t2, I528};
  auto task577 = make_shared<Task577>(tensor577, pindex);
  task576->add_dep(task577);
  task577->add_dep(task428);
  densityq->add_task(task577);

  vector<IndexRange> I529_index = {active_, virt_, closed_, active_};
  auto I529 = make_shared<Tensor>(I529_index);
  vector<shared_ptr<Tensor>> tensor578 = {I528, Gamma12_(), I529};
  auto task578 = make_shared<Task578>(tensor578, pindex);
  task577->add_dep(task578);
  task578->add_dep(task428);
  densityq->add_task(task578);

  vector<shared_ptr<Tensor>> tensor579 = {I529, t2};
  auto task579 = make_shared<Task579>(tensor579, pindex);
  task578->add_dep(task579);
  task579->add_dep(task428);
  densityq->add_task(task579);

  vector<IndexRange> I536_index = {virt_, active_};
  auto I536 = make_shared<Tensor>(I536_index);
  vector<shared_ptr<Tensor>> tensor580 = {den2, I536};
  auto task580 = make_shared<Task580>(tensor580, pindex);
  task580->add_dep(task428);
  densityq->add_task(task580);

  vector<IndexRange> I537_index = {active_, virt_};
  auto I537 = make_shared<Tensor>(I537_index);
  vector<shared_ptr<Tensor>> tensor581 = {I536, Gamma16_(), I537};
  auto task581 = make_shared<Task581>(tensor581, pindex);
  task580->add_dep(task581);
  task581->add_dep(task428);
  densityq->add_task(task581);

  vector<IndexRange> I538_index = {active_, closed_, virt_, closed_};
  auto I538 = make_shared<Tensor>(I538_index);
  vector<shared_ptr<Tensor>> tensor582 = {I537, t2, I538};
  auto task582 = make_shared<Task582>(tensor582, pindex);
  task581->add_dep(task582);
  task582->add_dep(task428);
  densityq->add_task(task582);

  vector<shared_ptr<Tensor>> tensor583 = {I538, t2};
  auto task583 = make_shared<Task583>(tensor583, pindex);
  task582->add_dep(task583);
  task583->add_dep(task428);
  densityq->add_task(task583);

  vector<IndexRange> I541_index = {active_, closed_, virt_, closed_};
  auto I541 = make_shared<Tensor>(I541_index);
  vector<shared_ptr<Tensor>> tensor584 = {I537, t2, I541};
  auto task584 = make_shared<Task584>(tensor584, pindex);
  task581->add_dep(task584);
  task584->add_dep(task428);
  densityq->add_task(task584);

  vector<shared_ptr<Tensor>> tensor585 = {I541, t2};
  auto task585 = make_shared<Task585>(tensor585, pindex);
  task584->add_dep(task585);
  task585->add_dep(task428);
  densityq->add_task(task585);

  vector<IndexRange> I542_index = {active_, virt_};
  auto I542 = make_shared<Tensor>(I542_index);
  vector<shared_ptr<Tensor>> tensor586 = {den2, I542};
  auto task586 = make_shared<Task586>(tensor586, pindex);
  task586->add_dep(task428);
  densityq->add_task(task586);

  vector<IndexRange> I543_index = {closed_, active_, active_, active_};
  auto I543 = make_shared<Tensor>(I543_index);
  vector<shared_ptr<Tensor>> tensor587 = {I542, t2, I543};
  auto task587 = make_shared<Task587>(tensor587, pindex);
  task586->add_dep(task587);
  task587->add_dep(task428);
  densityq->add_task(task587);

  vector<IndexRange> I544_index = {active_, active_, closed_, active_};
  auto I544 = make_shared<Tensor>(I544_index);
  vector<shared_ptr<Tensor>> tensor588 = {I543, Gamma28_(), I544};
  auto task588 = make_shared<Task588>(tensor588, pindex);
  task587->add_dep(task588);
  task588->add_dep(task428);
  densityq->add_task(task588);

  vector<shared_ptr<Tensor>> tensor589 = {I544, t2};
  auto task589 = make_shared<Task589>(tensor589, pindex);
  task588->add_dep(task589);
  task589->add_dep(task428);
  densityq->add_task(task589);

  vector<IndexRange> I545_index = {active_, closed_};
  auto I545 = make_shared<Tensor>(I545_index);
  vector<shared_ptr<Tensor>> tensor590 = {den2, I545};
  auto task590 = make_shared<Task590>(tensor590, pindex);
  task590->add_dep(task428);
  densityq->add_task(task590);

  vector<IndexRange> I546_index = {closed_, virt_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  vector<shared_ptr<Tensor>> tensor591 = {I545, t2, I546};
  auto task591 = make_shared<Task591>(tensor591, pindex);
  task590->add_dep(task591);
  task591->add_dep(task428);
  densityq->add_task(task591);

  vector<IndexRange> I547_index = {active_, closed_, virt_, active_};
  auto I547 = make_shared<Tensor>(I547_index);
  vector<shared_ptr<Tensor>> tensor592 = {I546, Gamma29_(), I547};
  auto task592 = make_shared<Task592>(tensor592, pindex);
  task591->add_dep(task592);
  task592->add_dep(task428);
  densityq->add_task(task592);

  vector<shared_ptr<Tensor>> tensor593 = {I547, t2};
  auto task593 = make_shared<Task593>(tensor593, pindex);
  task592->add_dep(task593);
  task593->add_dep(task428);
  densityq->add_task(task593);

  vector<IndexRange> I549_index = {closed_, virt_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  vector<shared_ptr<Tensor>> tensor594 = {I545, t2, I549};
  auto task594 = make_shared<Task594>(tensor594, pindex);
  task590->add_dep(task594);
  task594->add_dep(task428);
  densityq->add_task(task594);

  vector<IndexRange> I550_index = {active_, closed_, virt_, active_};
  auto I550 = make_shared<Tensor>(I550_index);
  vector<shared_ptr<Tensor>> tensor595 = {I549, Gamma7_(), I550};
  auto task595 = make_shared<Task595>(tensor595, pindex);
  task594->add_dep(task595);
  task595->add_dep(task428);
  densityq->add_task(task595);

  vector<shared_ptr<Tensor>> tensor596 = {I550, t2};
  auto task596 = make_shared<Task596>(tensor596, pindex);
  task595->add_dep(task596);
  task596->add_dep(task428);
  densityq->add_task(task596);

  vector<IndexRange> I588_index = {virt_, closed_, active_, active_};
  auto I588 = make_shared<Tensor>(I588_index);
  vector<shared_ptr<Tensor>> tensor597 = {I545, t2, I588};
  auto task597 = make_shared<Task597>(tensor597, pindex);
  task590->add_dep(task597);
  task597->add_dep(task428);
  densityq->add_task(task597);

  vector<IndexRange> I589_index = {active_, active_, virt_, closed_};
  auto I589 = make_shared<Tensor>(I589_index);
  vector<shared_ptr<Tensor>> tensor598 = {I588, Gamma7_(), I589};
  auto task598 = make_shared<Task598>(tensor598, pindex);
  task597->add_dep(task598);
  task598->add_dep(task428);
  densityq->add_task(task598);

  vector<shared_ptr<Tensor>> tensor599 = {I589, t2};
  auto task599 = make_shared<Task599>(tensor599, pindex);
  task598->add_dep(task599);
  task599->add_dep(task428);
  densityq->add_task(task599);

  vector<IndexRange> I591_index = {virt_, closed_, active_, active_};
  auto I591 = make_shared<Tensor>(I591_index);
  vector<shared_ptr<Tensor>> tensor600 = {I545, t2, I591};
  auto task600 = make_shared<Task600>(tensor600, pindex);
  task590->add_dep(task600);
  task600->add_dep(task428);
  densityq->add_task(task600);

  vector<IndexRange> I592_index = {active_, active_, virt_, closed_};
  auto I592 = make_shared<Tensor>(I592_index);
  vector<shared_ptr<Tensor>> tensor601 = {I591, Gamma7_(), I592};
  auto task601 = make_shared<Task601>(tensor601, pindex);
  task600->add_dep(task601);
  task601->add_dep(task428);
  densityq->add_task(task601);

  vector<shared_ptr<Tensor>> tensor602 = {I592, t2};
  auto task602 = make_shared<Task602>(tensor602, pindex);
  task601->add_dep(task602);
  task602->add_dep(task428);
  densityq->add_task(task602);

  vector<IndexRange> I740_index = {virt_, virt_, active_, active_};
  auto I740 = make_shared<Tensor>(I740_index);
  vector<shared_ptr<Tensor>> tensor603 = {I545, t2, I740};
  auto task603 = make_shared<Task603>(tensor603, pindex);
  task590->add_dep(task603);
  task603->add_dep(task428);
  densityq->add_task(task603);

  vector<IndexRange> I741_index = {virt_, active_, virt_, active_};
  auto I741 = make_shared<Tensor>(I741_index);
  vector<shared_ptr<Tensor>> tensor604 = {I740, Gamma60_(), I741};
  auto task604 = make_shared<Task604>(tensor604, pindex);
  task603->add_dep(task604);
  task604->add_dep(task428);
  densityq->add_task(task604);

  vector<shared_ptr<Tensor>> tensor605 = {I741, t2};
  auto task605 = make_shared<Task605>(tensor605, pindex);
  task604->add_dep(task605);
  task605->add_dep(task428);
  densityq->add_task(task605);

  vector<IndexRange> I557_index = {virt_, virt_};
  auto I557 = make_shared<Tensor>(I557_index);
  vector<shared_ptr<Tensor>> tensor606 = {den2, I557};
  auto task606 = make_shared<Task606>(tensor606, pindex);
  task606->add_dep(task428);
  densityq->add_task(task606);

  vector<IndexRange> I558_index = {closed_, virt_, active_, active_};
  auto I558 = make_shared<Tensor>(I558_index);
  vector<shared_ptr<Tensor>> tensor607 = {I557, t2, I558};
  auto task607 = make_shared<Task607>(tensor607, pindex);
  task606->add_dep(task607);
  task607->add_dep(task428);
  densityq->add_task(task607);

  vector<IndexRange> I559_index = {active_, closed_, virt_, active_};
  auto I559 = make_shared<Tensor>(I559_index);
  vector<shared_ptr<Tensor>> tensor608 = {I558, Gamma32_(), I559};
  auto task608 = make_shared<Task608>(tensor608, pindex);
  task607->add_dep(task608);
  task608->add_dep(task428);
  densityq->add_task(task608);

  vector<shared_ptr<Tensor>> tensor609 = {I559, t2};
  auto task609 = make_shared<Task609>(tensor609, pindex);
  task608->add_dep(task609);
  task609->add_dep(task428);
  densityq->add_task(task609);

  vector<IndexRange> I567_index = {closed_, virt_, active_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  vector<shared_ptr<Tensor>> tensor610 = {I557, t2, I567};
  auto task610 = make_shared<Task610>(tensor610, pindex);
  task606->add_dep(task610);
  task610->add_dep(task428);
  densityq->add_task(task610);

  vector<IndexRange> I568_index = {active_, closed_, virt_, active_};
  auto I568 = make_shared<Tensor>(I568_index);
  vector<shared_ptr<Tensor>> tensor611 = {I567, Gamma35_(), I568};
  auto task611 = make_shared<Task611>(tensor611, pindex);
  task610->add_dep(task611);
  task611->add_dep(task428);
  densityq->add_task(task611);

  vector<shared_ptr<Tensor>> tensor612 = {I568, t2};
  auto task612 = make_shared<Task612>(tensor612, pindex);
  task611->add_dep(task612);
  task612->add_dep(task428);
  densityq->add_task(task612);

  vector<IndexRange> I572_index = {virt_, closed_};
  auto I572 = make_shared<Tensor>(I572_index);
  vector<shared_ptr<Tensor>> tensor613 = {den2, I572};
  auto task613 = make_shared<Task613>(tensor613, pindex);
  task613->add_dep(task428);
  densityq->add_task(task613);

  vector<IndexRange> I573_index = {closed_, virt_};
  auto I573 = make_shared<Tensor>(I573_index);
  vector<shared_ptr<Tensor>> tensor614 = {I572, t2, I573};
  auto task614 = make_shared<Task614>(tensor614, pindex);
  task613->add_dep(task614);
  task614->add_dep(task428);
  densityq->add_task(task614);

  vector<IndexRange> I574_index = {active_, closed_, virt_, active_};
  auto I574 = make_shared<Tensor>(I574_index);
  vector<shared_ptr<Tensor>> tensor615 = {I573, Gamma38_(), I574};
  auto task615 = make_shared<Task615>(tensor615, pindex);
  task614->add_dep(task615);
  task615->add_dep(task428);
  densityq->add_task(task615);

  vector<shared_ptr<Tensor>> tensor616 = {I574, t2};
  auto task616 = make_shared<Task616>(tensor616, pindex);
  task615->add_dep(task616);
  task616->add_dep(task428);
  densityq->add_task(task616);

  vector<IndexRange> I576_index = {closed_, virt_};
  auto I576 = make_shared<Tensor>(I576_index);
  vector<shared_ptr<Tensor>> tensor617 = {I572, t2, I576};
  auto task617 = make_shared<Task617>(tensor617, pindex);
  task613->add_dep(task617);
  task617->add_dep(task428);
  densityq->add_task(task617);

  vector<IndexRange> I577_index = {active_, closed_, virt_, active_};
  auto I577 = make_shared<Tensor>(I577_index);
  vector<shared_ptr<Tensor>> tensor618 = {I576, Gamma38_(), I577};
  auto task618 = make_shared<Task618>(tensor618, pindex);
  task617->add_dep(task618);
  task618->add_dep(task428);
  densityq->add_task(task618);

  vector<shared_ptr<Tensor>> tensor619 = {I577, t2};
  auto task619 = make_shared<Task619>(tensor619, pindex);
  task618->add_dep(task619);
  task619->add_dep(task428);
  densityq->add_task(task619);

  vector<IndexRange> I615_index = {virt_, closed_};
  auto I615 = make_shared<Tensor>(I615_index);
  vector<shared_ptr<Tensor>> tensor620 = {I572, t2, I615};
  auto task620 = make_shared<Task620>(tensor620, pindex);
  task613->add_dep(task620);
  task620->add_dep(task428);
  densityq->add_task(task620);

  vector<IndexRange> I616_index = {active_, active_, virt_, closed_};
  auto I616 = make_shared<Tensor>(I616_index);
  vector<shared_ptr<Tensor>> tensor621 = {I615, Gamma38_(), I616};
  auto task621 = make_shared<Task621>(tensor621, pindex);
  task620->add_dep(task621);
  task621->add_dep(task428);
  densityq->add_task(task621);

  vector<shared_ptr<Tensor>> tensor622 = {I616, t2};
  auto task622 = make_shared<Task622>(tensor622, pindex);
  task621->add_dep(task622);
  task622->add_dep(task428);
  densityq->add_task(task622);

  vector<IndexRange> I618_index = {virt_, closed_};
  auto I618 = make_shared<Tensor>(I618_index);
  vector<shared_ptr<Tensor>> tensor623 = {I572, t2, I618};
  auto task623 = make_shared<Task623>(tensor623, pindex);
  task613->add_dep(task623);
  task623->add_dep(task428);
  densityq->add_task(task623);

  vector<IndexRange> I619_index = {active_, active_, virt_, closed_};
  auto I619 = make_shared<Tensor>(I619_index);
  vector<shared_ptr<Tensor>> tensor624 = {I618, Gamma38_(), I619};
  auto task624 = make_shared<Task624>(tensor624, pindex);
  task623->add_dep(task624);
  task624->add_dep(task428);
  densityq->add_task(task624);

  vector<shared_ptr<Tensor>> tensor625 = {I619, t2};
  auto task625 = make_shared<Task625>(tensor625, pindex);
  task624->add_dep(task625);
  task625->add_dep(task428);
  densityq->add_task(task625);

  vector<IndexRange> I584_index = {active_, virt_};
  auto I584 = make_shared<Tensor>(I584_index);
  vector<shared_ptr<Tensor>> tensor626 = {den2, I584};
  auto task626 = make_shared<Task626>(tensor626, pindex);
  task626->add_dep(task428);
  densityq->add_task(task626);

  vector<IndexRange> I585_index = {closed_, active_, active_, active_};
  auto I585 = make_shared<Tensor>(I585_index);
  vector<shared_ptr<Tensor>> tensor627 = {I584, t2, I585};
  auto task627 = make_shared<Task627>(tensor627, pindex);
  task626->add_dep(task627);
  task627->add_dep(task428);
  densityq->add_task(task627);

  vector<IndexRange> I586_index = {active_, active_, closed_, active_};
  auto I586 = make_shared<Tensor>(I586_index);
  vector<shared_ptr<Tensor>> tensor628 = {I585, Gamma6_(), I586};
  auto task628 = make_shared<Task628>(tensor628, pindex);
  task627->add_dep(task628);
  task628->add_dep(task428);
  densityq->add_task(task628);

  vector<shared_ptr<Tensor>> tensor629 = {I586, t2};
  auto task629 = make_shared<Task629>(tensor629, pindex);
  task628->add_dep(task629);
  task629->add_dep(task428);
  densityq->add_task(task629);

  vector<IndexRange> I737_index = {virt_, active_, active_, active_};
  auto I737 = make_shared<Tensor>(I737_index);
  vector<shared_ptr<Tensor>> tensor630 = {I584, t2, I737};
  auto task630 = make_shared<Task630>(tensor630, pindex);
  task626->add_dep(task630);
  task630->add_dep(task428);
  densityq->add_task(task630);

  vector<IndexRange> I738_index = {active_, virt_, active_, active_};
  auto I738 = make_shared<Tensor>(I738_index);
  vector<shared_ptr<Tensor>> tensor631 = {I737, Gamma59_(), I738};
  auto task631 = make_shared<Task631>(tensor631, pindex);
  task630->add_dep(task631);
  task631->add_dep(task428);
  densityq->add_task(task631);

  vector<shared_ptr<Tensor>> tensor632 = {I738, t2};
  auto task632 = make_shared<Task632>(tensor632, pindex);
  task631->add_dep(task632);
  task632->add_dep(task428);
  densityq->add_task(task632);

  vector<IndexRange> I596_index = {closed_, closed_};
  auto I596 = make_shared<Tensor>(I596_index);
  vector<shared_ptr<Tensor>> tensor633 = {den2, I596};
  auto task633 = make_shared<Task633>(tensor633, pindex);
  task633->add_dep(task428);
  densityq->add_task(task633);

  vector<IndexRange> I597_index = {virt_, closed_, active_, active_};
  auto I597 = make_shared<Tensor>(I597_index);
  vector<shared_ptr<Tensor>> tensor634 = {I596, t2, I597};
  auto task634 = make_shared<Task634>(tensor634, pindex);
  task633->add_dep(task634);
  task634->add_dep(task428);
  densityq->add_task(task634);

  vector<IndexRange> I598_index = {active_, active_, virt_, closed_};
  auto I598 = make_shared<Tensor>(I598_index);
  vector<shared_ptr<Tensor>> tensor635 = {I597, Gamma35_(), I598};
  auto task635 = make_shared<Task635>(tensor635, pindex);
  task634->add_dep(task635);
  task635->add_dep(task428);
  densityq->add_task(task635);

  vector<shared_ptr<Tensor>> tensor636 = {I598, t2};
  auto task636 = make_shared<Task636>(tensor636, pindex);
  task635->add_dep(task636);
  task636->add_dep(task428);
  densityq->add_task(task636);

  vector<IndexRange> I606_index = {virt_, closed_, active_, active_};
  auto I606 = make_shared<Tensor>(I606_index);
  vector<shared_ptr<Tensor>> tensor637 = {I596, t2, I606};
  auto task637 = make_shared<Task637>(tensor637, pindex);
  task633->add_dep(task637);
  task637->add_dep(task428);
  densityq->add_task(task637);

  vector<IndexRange> I607_index = {active_, active_, virt_, closed_};
  auto I607 = make_shared<Tensor>(I607_index);
  vector<shared_ptr<Tensor>> tensor638 = {I606, Gamma35_(), I607};
  auto task638 = make_shared<Task638>(tensor638, pindex);
  task637->add_dep(task638);
  task638->add_dep(task428);
  densityq->add_task(task638);

  vector<shared_ptr<Tensor>> tensor639 = {I607, t2};
  auto task639 = make_shared<Task639>(tensor639, pindex);
  task638->add_dep(task639);
  task639->add_dep(task428);
  densityq->add_task(task639);

  vector<IndexRange> I599_index = {virt_, virt_};
  auto I599 = make_shared<Tensor>(I599_index);
  vector<shared_ptr<Tensor>> tensor640 = {den2, I599};
  auto task640 = make_shared<Task640>(tensor640, pindex);
  task640->add_dep(task428);
  densityq->add_task(task640);

  vector<IndexRange> I600_index = {virt_, closed_, active_, active_};
  auto I600 = make_shared<Tensor>(I600_index);
  vector<shared_ptr<Tensor>> tensor641 = {I599, t2, I600};
  auto task641 = make_shared<Task641>(tensor641, pindex);
  task640->add_dep(task641);
  task641->add_dep(task428);
  densityq->add_task(task641);

  vector<IndexRange> I601_index = {active_, active_, virt_, closed_};
  auto I601 = make_shared<Tensor>(I601_index);
  vector<shared_ptr<Tensor>> tensor642 = {I600, Gamma35_(), I601};
  auto task642 = make_shared<Task642>(tensor642, pindex);
  task641->add_dep(task642);
  task642->add_dep(task428);
  densityq->add_task(task642);

  vector<shared_ptr<Tensor>> tensor643 = {I601, t2};
  auto task643 = make_shared<Task643>(tensor643, pindex);
  task642->add_dep(task643);
  task643->add_dep(task428);
  densityq->add_task(task643);

  vector<IndexRange> I609_index = {virt_, closed_, active_, active_};
  auto I609 = make_shared<Tensor>(I609_index);
  vector<shared_ptr<Tensor>> tensor644 = {I599, t2, I609};
  auto task644 = make_shared<Task644>(tensor644, pindex);
  task640->add_dep(task644);
  task644->add_dep(task428);
  densityq->add_task(task644);

  vector<IndexRange> I610_index = {active_, active_, virt_, closed_};
  auto I610 = make_shared<Tensor>(I610_index);
  vector<shared_ptr<Tensor>> tensor645 = {I609, Gamma35_(), I610};
  auto task645 = make_shared<Task645>(tensor645, pindex);
  task644->add_dep(task645);
  task645->add_dep(task428);
  densityq->add_task(task645);

  vector<shared_ptr<Tensor>> tensor646 = {I610, t2};
  auto task646 = make_shared<Task646>(tensor646, pindex);
  task645->add_dep(task646);
  task646->add_dep(task428);
  densityq->add_task(task646);

  vector<IndexRange> I746_index = {virt_, virt_, active_, active_};
  auto I746 = make_shared<Tensor>(I746_index);
  vector<shared_ptr<Tensor>> tensor647 = {I599, t2, I746};
  auto task647 = make_shared<Task647>(tensor647, pindex);
  task640->add_dep(task647);
  task647->add_dep(task428);
  densityq->add_task(task647);

  vector<IndexRange> I747_index = {virt_, active_, virt_, active_};
  auto I747 = make_shared<Tensor>(I747_index);
  vector<shared_ptr<Tensor>> tensor648 = {I746, Gamma60_(), I747};
  auto task648 = make_shared<Task648>(tensor648, pindex);
  task647->add_dep(task648);
  task648->add_dep(task428);
  densityq->add_task(task648);

  vector<shared_ptr<Tensor>> tensor649 = {I747, t2};
  auto task649 = make_shared<Task649>(tensor649, pindex);
  task648->add_dep(task649);
  task649->add_dep(task428);
  densityq->add_task(task649);

  vector<IndexRange> I611_index = {active_, closed_};
  auto I611 = make_shared<Tensor>(I611_index);
  vector<shared_ptr<Tensor>> tensor650 = {den2, I611};
  auto task650 = make_shared<Task650>(tensor650, pindex);
  task650->add_dep(task428);
  densityq->add_task(task650);

  vector<IndexRange> I612_index = {virt_, active_, active_, active_};
  auto I612 = make_shared<Tensor>(I612_index);
  vector<shared_ptr<Tensor>> tensor651 = {I611, t2, I612};
  auto task651 = make_shared<Task651>(tensor651, pindex);
  task650->add_dep(task651);
  task651->add_dep(task428);
  densityq->add_task(task651);

  vector<IndexRange> I613_index = {active_, virt_, active_, active_};
  auto I613 = make_shared<Tensor>(I613_index);
  vector<shared_ptr<Tensor>> tensor652 = {I612, Gamma51_(), I613};
  auto task652 = make_shared<Task652>(tensor652, pindex);
  task651->add_dep(task652);
  task652->add_dep(task428);
  densityq->add_task(task652);

  vector<shared_ptr<Tensor>> tensor653 = {I613, t2};
  auto task653 = make_shared<Task653>(tensor653, pindex);
  task652->add_dep(task653);
  task653->add_dep(task428);
  densityq->add_task(task653);

  vector<IndexRange> I635_index = {virt_, virt_};
  auto I635 = make_shared<Tensor>(I635_index);
  vector<shared_ptr<Tensor>> tensor654 = {den2, I635};
  auto task654 = make_shared<Task654>(tensor654, pindex);
  task654->add_dep(task428);
  densityq->add_task(task654);

  vector<IndexRange> I636_index = {virt_, active_, active_, active_};
  auto I636 = make_shared<Tensor>(I636_index);
  vector<shared_ptr<Tensor>> tensor655 = {I635, t2, I636};
  auto task655 = make_shared<Task655>(tensor655, pindex);
  task654->add_dep(task655);
  task655->add_dep(task428);
  densityq->add_task(task655);

  vector<IndexRange> I637_index = {active_, active_, virt_, active_};
  auto I637 = make_shared<Tensor>(I637_index);
  vector<shared_ptr<Tensor>> tensor656 = {I636, Gamma59_(), I637};
  auto task656 = make_shared<Task656>(tensor656, pindex);
  task655->add_dep(task656);
  task656->add_dep(task428);
  densityq->add_task(task656);

  vector<shared_ptr<Tensor>> tensor657 = {I637, t2};
  auto task657 = make_shared<Task657>(tensor657, pindex);
  task656->add_dep(task657);
  task657->add_dep(task428);
  densityq->add_task(task657);

  vector<IndexRange> I647_index = {virt_, active_};
  auto I647 = make_shared<Tensor>(I647_index);
  vector<shared_ptr<Tensor>> tensor658 = {den2, I647};
  auto task658 = make_shared<Task658>(tensor658, pindex);
  task658->add_dep(task428);
  densityq->add_task(task658);

  vector<IndexRange> I648_index = {virt_, active_};
  auto I648 = make_shared<Tensor>(I648_index);
  vector<shared_ptr<Tensor>> tensor659 = {I647, Gamma16_(), I648};
  auto task659 = make_shared<Task659>(tensor659, pindex);
  task658->add_dep(task659);
  task659->add_dep(task428);
  densityq->add_task(task659);

  vector<IndexRange> I649_index = {virt_, closed_, virt_, closed_};
  auto I649 = make_shared<Tensor>(I649_index);
  vector<shared_ptr<Tensor>> tensor660 = {I648, t2, I649};
  auto task660 = make_shared<Task660>(tensor660, pindex);
  task659->add_dep(task660);
  task660->add_dep(task428);
  densityq->add_task(task660);

  vector<shared_ptr<Tensor>> tensor661 = {I649, t2};
  auto task661 = make_shared<Task661>(tensor661, pindex);
  task660->add_dep(task661);
  task661->add_dep(task428);
  densityq->add_task(task661);

  vector<IndexRange> I650_index = {virt_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  vector<shared_ptr<Tensor>> tensor662 = {den2, I650};
  auto task662 = make_shared<Task662>(tensor662, pindex);
  task662->add_dep(task428);
  densityq->add_task(task662);

  vector<IndexRange> I651_index = {virt_, active_};
  auto I651 = make_shared<Tensor>(I651_index);
  vector<shared_ptr<Tensor>> tensor663 = {I650, Gamma16_(), I651};
  auto task663 = make_shared<Task663>(tensor663, pindex);
  task662->add_dep(task663);
  task663->add_dep(task428);
  densityq->add_task(task663);

  vector<IndexRange> I652_index = {virt_, closed_, virt_, closed_};
  auto I652 = make_shared<Tensor>(I652_index);
  vector<shared_ptr<Tensor>> tensor664 = {I651, t2, I652};
  auto task664 = make_shared<Task664>(tensor664, pindex);
  task663->add_dep(task664);
  task664->add_dep(task428);
  densityq->add_task(task664);

  vector<shared_ptr<Tensor>> tensor665 = {I652, t2};
  auto task665 = make_shared<Task665>(tensor665, pindex);
  task664->add_dep(task665);
  task665->add_dep(task428);
  densityq->add_task(task665);

  vector<IndexRange> I656_index = {virt_, closed_};
  auto I656 = make_shared<Tensor>(I656_index);
  vector<shared_ptr<Tensor>> tensor666 = {den2, I656};
  auto task666 = make_shared<Task666>(tensor666, pindex);
  task666->add_dep(task428);
  densityq->add_task(task666);

  vector<IndexRange> I657_index = {virt_, closed_};
  auto I657 = make_shared<Tensor>(I657_index);
  vector<shared_ptr<Tensor>> tensor667 = {I656, t2, I657};
  auto task667 = make_shared<Task667>(tensor667, pindex);
  task666->add_dep(task667);
  task667->add_dep(task428);
  densityq->add_task(task667);

  vector<IndexRange> I658_index = {active_, virt_, closed_, active_};
  auto I658 = make_shared<Tensor>(I658_index);
  vector<shared_ptr<Tensor>> tensor668 = {I657, Gamma38_(), I658};
  auto task668 = make_shared<Task668>(tensor668, pindex);
  task667->add_dep(task668);
  task668->add_dep(task428);
  densityq->add_task(task668);

  vector<shared_ptr<Tensor>> tensor669 = {I658, t2};
  auto task669 = make_shared<Task669>(tensor669, pindex);
  task668->add_dep(task669);
  task669->add_dep(task428);
  densityq->add_task(task669);

  vector<IndexRange> I665_index = {active_, active_};
  auto I665 = make_shared<Tensor>(I665_index);
  vector<shared_ptr<Tensor>> tensor670 = {den2, I665};
  auto task670 = make_shared<Task670>(tensor670, pindex);
  task670->add_dep(task428);
  densityq->add_task(task670);

  vector<IndexRange> I666_index;
  auto I666 = make_shared<Tensor>(I666_index);
  vector<shared_ptr<Tensor>> tensor671 = {I665, Gamma38_(), I666};
  auto task671 = make_shared<Task671>(tensor671, pindex);
  task670->add_dep(task671);
  task671->add_dep(task428);
  densityq->add_task(task671);

  vector<IndexRange> I667_index = {virt_, closed_, virt_, closed_};
  auto I667 = make_shared<Tensor>(I667_index);
  vector<shared_ptr<Tensor>> tensor672 = {I666, t2, I667};
  auto task672 = make_shared<Task672>(tensor672, pindex);
  task671->add_dep(task672);
  task672->add_dep(task428);
  densityq->add_task(task672);

  vector<shared_ptr<Tensor>> tensor673 = {I667, t2};
  auto task673 = make_shared<Task673>(tensor673, pindex);
  task672->add_dep(task673);
  task673->add_dep(task428);
  densityq->add_task(task673);

  vector<IndexRange> I670_index = {virt_, closed_, virt_, closed_};
  auto I670 = make_shared<Tensor>(I670_index);
  vector<shared_ptr<Tensor>> tensor674 = {I666, t2, I670};
  auto task674 = make_shared<Task674>(tensor674, pindex);
  task671->add_dep(task674);
  task674->add_dep(task428);
  densityq->add_task(task674);

  vector<shared_ptr<Tensor>> tensor675 = {I670, t2};
  auto task675 = make_shared<Task675>(tensor675, pindex);
  task674->add_dep(task675);
  task675->add_dep(task428);
  densityq->add_task(task675);

  vector<IndexRange> I671_index = {closed_, closed_};
  auto I671 = make_shared<Tensor>(I671_index);
  vector<shared_ptr<Tensor>> tensor676 = {den2, I671};
  auto task676 = make_shared<Task676>(tensor676, pindex);
  task676->add_dep(task428);
  densityq->add_task(task676);

  vector<IndexRange> I672_index = {virt_, closed_, virt_, closed_};
  auto I672 = make_shared<Tensor>(I672_index);
  vector<shared_ptr<Tensor>> tensor677 = {I671, t2, I672};
  auto task677 = make_shared<Task677>(tensor677, pindex);
  task676->add_dep(task677);
  task677->add_dep(task428);
  densityq->add_task(task677);

  vector<shared_ptr<Tensor>> tensor678 = {I672, t2};
  auto task678 = make_shared<Task678>(tensor678, pindex);
  task677->add_dep(task678);
  task678->add_dep(task428);
  densityq->add_task(task678);

  vector<IndexRange> I674_index = {virt_, closed_, virt_, closed_};
  auto I674 = make_shared<Tensor>(I674_index);
  vector<shared_ptr<Tensor>> tensor679 = {I671, t2, I674};
  auto task679 = make_shared<Task679>(tensor679, pindex);
  task676->add_dep(task679);
  task679->add_dep(task428);
  densityq->add_task(task679);

  vector<shared_ptr<Tensor>> tensor680 = {I674, t2};
  auto task680 = make_shared<Task680>(tensor680, pindex);
  task679->add_dep(task680);
  task680->add_dep(task428);
  densityq->add_task(task680);

  vector<IndexRange> I675_index = {virt_, virt_};
  auto I675 = make_shared<Tensor>(I675_index);
  vector<shared_ptr<Tensor>> tensor681 = {den2, I675};
  auto task681 = make_shared<Task681>(tensor681, pindex);
  task681->add_dep(task428);
  densityq->add_task(task681);

  vector<IndexRange> I676_index = {virt_, closed_, virt_, closed_};
  auto I676 = make_shared<Tensor>(I676_index);
  vector<shared_ptr<Tensor>> tensor682 = {I675, t2, I676};
  auto task682 = make_shared<Task682>(tensor682, pindex);
  task681->add_dep(task682);
  task682->add_dep(task428);
  densityq->add_task(task682);

  vector<shared_ptr<Tensor>> tensor683 = {I676, t2};
  auto task683 = make_shared<Task683>(tensor683, pindex);
  task682->add_dep(task683);
  task683->add_dep(task428);
  densityq->add_task(task683);

  vector<IndexRange> I678_index = {virt_, closed_, virt_, closed_};
  auto I678 = make_shared<Tensor>(I678_index);
  vector<shared_ptr<Tensor>> tensor684 = {I675, t2, I678};
  auto task684 = make_shared<Task684>(tensor684, pindex);
  task681->add_dep(task684);
  task684->add_dep(task428);
  densityq->add_task(task684);

  vector<shared_ptr<Tensor>> tensor685 = {I678, t2};
  auto task685 = make_shared<Task685>(tensor685, pindex);
  task684->add_dep(task685);
  task685->add_dep(task428);
  densityq->add_task(task685);

  vector<IndexRange> I679_index = {closed_, active_};
  auto I679 = make_shared<Tensor>(I679_index);
  vector<shared_ptr<Tensor>> tensor686 = {den2, I679};
  auto task686 = make_shared<Task686>(tensor686, pindex);
  task686->add_dep(task428);
  densityq->add_task(task686);

  vector<IndexRange> I680_index = {closed_, active_};
  auto I680 = make_shared<Tensor>(I680_index);
  vector<shared_ptr<Tensor>> tensor687 = {I679, Gamma38_(), I680};
  auto task687 = make_shared<Task687>(tensor687, pindex);
  task686->add_dep(task687);
  task687->add_dep(task428);
  densityq->add_task(task687);

  vector<IndexRange> I681_index = {virt_, closed_, virt_, closed_};
  auto I681 = make_shared<Tensor>(I681_index);
  vector<shared_ptr<Tensor>> tensor688 = {I680, t2, I681};
  auto task688 = make_shared<Task688>(tensor688, pindex);
  task687->add_dep(task688);
  task688->add_dep(task428);
  densityq->add_task(task688);

  vector<shared_ptr<Tensor>> tensor689 = {I681, t2};
  auto task689 = make_shared<Task689>(tensor689, pindex);
  task688->add_dep(task689);
  task689->add_dep(task428);
  densityq->add_task(task689);

  vector<IndexRange> I684_index = {virt_, closed_, virt_, closed_};
  auto I684 = make_shared<Tensor>(I684_index);
  vector<shared_ptr<Tensor>> tensor690 = {I680, t2, I684};
  auto task690 = make_shared<Task690>(tensor690, pindex);
  task687->add_dep(task690);
  task690->add_dep(task428);
  densityq->add_task(task690);

  vector<shared_ptr<Tensor>> tensor691 = {I684, t2};
  auto task691 = make_shared<Task691>(tensor691, pindex);
  task690->add_dep(task691);
  task691->add_dep(task428);
  densityq->add_task(task691);

  vector<IndexRange> I685_index = {active_, virt_};
  auto I685 = make_shared<Tensor>(I685_index);
  vector<shared_ptr<Tensor>> tensor692 = {den2, I685};
  auto task692 = make_shared<Task692>(tensor692, pindex);
  task692->add_dep(task428);
  densityq->add_task(task692);

  vector<IndexRange> I686_index = {virt_, closed_, active_, active_};
  auto I686 = make_shared<Tensor>(I686_index);
  vector<shared_ptr<Tensor>> tensor693 = {I685, t2, I686};
  auto task693 = make_shared<Task693>(tensor693, pindex);
  task692->add_dep(task693);
  task693->add_dep(task428);
  densityq->add_task(task693);

  vector<IndexRange> I687_index = {active_, virt_, closed_, active_};
  auto I687 = make_shared<Tensor>(I687_index);
  vector<shared_ptr<Tensor>> tensor694 = {I686, Gamma35_(), I687};
  auto task694 = make_shared<Task694>(tensor694, pindex);
  task693->add_dep(task694);
  task694->add_dep(task428);
  densityq->add_task(task694);

  vector<shared_ptr<Tensor>> tensor695 = {I687, t2};
  auto task695 = make_shared<Task695>(tensor695, pindex);
  task694->add_dep(task695);
  task695->add_dep(task428);
  densityq->add_task(task695);

  vector<IndexRange> I688_index = {active_, virt_};
  auto I688 = make_shared<Tensor>(I688_index);
  vector<shared_ptr<Tensor>> tensor696 = {den2, I688};
  auto task696 = make_shared<Task696>(tensor696, pindex);
  task696->add_dep(task428);
  densityq->add_task(task696);

  vector<IndexRange> I689_index = {virt_, closed_, active_, active_};
  auto I689 = make_shared<Tensor>(I689_index);
  vector<shared_ptr<Tensor>> tensor697 = {I688, t2, I689};
  auto task697 = make_shared<Task697>(tensor697, pindex);
  task696->add_dep(task697);
  task697->add_dep(task428);
  densityq->add_task(task697);

  vector<IndexRange> I690_index = {active_, virt_, closed_, active_};
  auto I690 = make_shared<Tensor>(I690_index);
  vector<shared_ptr<Tensor>> tensor698 = {I689, Gamma32_(), I690};
  auto task698 = make_shared<Task698>(tensor698, pindex);
  task697->add_dep(task698);
  task698->add_dep(task428);
  densityq->add_task(task698);

  vector<shared_ptr<Tensor>> tensor699 = {I690, t2};
  auto task699 = make_shared<Task699>(tensor699, pindex);
  task698->add_dep(task699);
  task699->add_dep(task428);
  densityq->add_task(task699);

  vector<IndexRange> I696_index = {closed_, virt_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  vector<shared_ptr<Tensor>> tensor700 = {I689, Gamma35_(), I696};
  auto task700 = make_shared<Task700>(tensor700, pindex);
  task697->add_dep(task700);
  task700->add_dep(task428);
  densityq->add_task(task700);

  vector<shared_ptr<Tensor>> tensor701 = {I696, t2};
  auto task701 = make_shared<Task701>(tensor701, pindex);
  task700->add_dep(task701);
  task701->add_dep(task428);
  densityq->add_task(task701);

  vector<IndexRange> I697_index = {closed_, virt_};
  auto I697 = make_shared<Tensor>(I697_index);
  vector<shared_ptr<Tensor>> tensor702 = {den2, I697};
  auto task702 = make_shared<Task702>(tensor702, pindex);
  task702->add_dep(task428);
  densityq->add_task(task702);

  vector<IndexRange> I698_index = {virt_, active_};
  auto I698 = make_shared<Tensor>(I698_index);
  vector<shared_ptr<Tensor>> tensor703 = {I697, t2, I698};
  auto task703 = make_shared<Task703>(tensor703, pindex);
  task702->add_dep(task703);
  task703->add_dep(task428);
  densityq->add_task(task703);

  vector<IndexRange> I699_index = {active_, virt_, active_, active_};
  auto I699 = make_shared<Tensor>(I699_index);
  vector<shared_ptr<Tensor>> tensor704 = {I698, Gamma60_(), I699};
  auto task704 = make_shared<Task704>(tensor704, pindex);
  task703->add_dep(task704);
  task704->add_dep(task428);
  densityq->add_task(task704);

  vector<shared_ptr<Tensor>> tensor705 = {I699, t2};
  auto task705 = make_shared<Task705>(tensor705, pindex);
  task704->add_dep(task705);
  task705->add_dep(task428);
  densityq->add_task(task705);

  vector<IndexRange> I700_index = {virt_, closed_};
  auto I700 = make_shared<Tensor>(I700_index);
  vector<shared_ptr<Tensor>> tensor706 = {den2, I700};
  auto task706 = make_shared<Task706>(tensor706, pindex);
  task706->add_dep(task428);
  densityq->add_task(task706);

  vector<IndexRange> I701_index = {virt_, active_};
  auto I701 = make_shared<Tensor>(I701_index);
  vector<shared_ptr<Tensor>> tensor707 = {I700, t2, I701};
  auto task707 = make_shared<Task707>(tensor707, pindex);
  task706->add_dep(task707);
  task707->add_dep(task428);
  densityq->add_task(task707);

  vector<IndexRange> I702_index = {active_, virt_, active_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  vector<shared_ptr<Tensor>> tensor708 = {I701, Gamma60_(), I702};
  auto task708 = make_shared<Task708>(tensor708, pindex);
  task707->add_dep(task708);
  task708->add_dep(task428);
  densityq->add_task(task708);

  vector<shared_ptr<Tensor>> tensor709 = {I702, t2};
  auto task709 = make_shared<Task709>(tensor709, pindex);
  task708->add_dep(task709);
  task709->add_dep(task428);
  densityq->add_task(task709);

  vector<IndexRange> I703_index = {closed_, active_};
  auto I703 = make_shared<Tensor>(I703_index);
  vector<shared_ptr<Tensor>> tensor710 = {den2, I703};
  auto task710 = make_shared<Task710>(tensor710, pindex);
  task710->add_dep(task428);
  densityq->add_task(task710);

  vector<IndexRange> I704_index = {active_, closed_};
  auto I704 = make_shared<Tensor>(I704_index);
  vector<shared_ptr<Tensor>> tensor711 = {I703, Gamma38_(), I704};
  auto task711 = make_shared<Task711>(tensor711, pindex);
  task710->add_dep(task711);
  task711->add_dep(task428);
  densityq->add_task(task711);

  vector<IndexRange> I705_index = {virt_, closed_, virt_, active_};
  auto I705 = make_shared<Tensor>(I705_index);
  vector<shared_ptr<Tensor>> tensor712 = {I704, t2, I705};
  auto task712 = make_shared<Task712>(tensor712, pindex);
  task711->add_dep(task712);
  task712->add_dep(task428);
  densityq->add_task(task712);

  vector<shared_ptr<Tensor>> tensor713 = {I705, t2};
  auto task713 = make_shared<Task713>(tensor713, pindex);
  task712->add_dep(task713);
  task713->add_dep(task428);
  densityq->add_task(task713);

  vector<IndexRange> I708_index = {virt_, closed_, virt_, active_};
  auto I708 = make_shared<Tensor>(I708_index);
  vector<shared_ptr<Tensor>> tensor714 = {I704, t2, I708};
  auto task714 = make_shared<Task714>(tensor714, pindex);
  task711->add_dep(task714);
  task714->add_dep(task428);
  densityq->add_task(task714);

  vector<shared_ptr<Tensor>> tensor715 = {I708, t2};
  auto task715 = make_shared<Task715>(tensor715, pindex);
  task714->add_dep(task715);
  task715->add_dep(task428);
  densityq->add_task(task715);

  vector<IndexRange> I715_index = {closed_, closed_};
  auto I715 = make_shared<Tensor>(I715_index);
  vector<shared_ptr<Tensor>> tensor716 = {den2, I715};
  auto task716 = make_shared<Task716>(tensor716, pindex);
  task716->add_dep(task428);
  densityq->add_task(task716);

  vector<IndexRange> I716_index = {virt_, closed_, virt_, active_};
  auto I716 = make_shared<Tensor>(I716_index);
  vector<shared_ptr<Tensor>> tensor717 = {I715, t2, I716};
  auto task717 = make_shared<Task717>(tensor717, pindex);
  task716->add_dep(task717);
  task717->add_dep(task428);
  densityq->add_task(task717);

  vector<IndexRange> I717_index = {virt_, closed_, virt_, active_};
  auto I717 = make_shared<Tensor>(I717_index);
  vector<shared_ptr<Tensor>> tensor718 = {I716, Gamma38_(), I717};
  auto task718 = make_shared<Task718>(tensor718, pindex);
  task717->add_dep(task718);
  task718->add_dep(task428);
  densityq->add_task(task718);

  vector<shared_ptr<Tensor>> tensor719 = {I717, t2};
  auto task719 = make_shared<Task719>(tensor719, pindex);
  task718->add_dep(task719);
  task719->add_dep(task428);
  densityq->add_task(task719);

  vector<IndexRange> I719_index = {virt_, closed_, virt_, active_};
  auto I719 = make_shared<Tensor>(I719_index);
  vector<shared_ptr<Tensor>> tensor720 = {I715, t2, I719};
  auto task720 = make_shared<Task720>(tensor720, pindex);
  task716->add_dep(task720);
  task720->add_dep(task428);
  densityq->add_task(task720);

  vector<IndexRange> I720_index = {virt_, closed_, virt_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  vector<shared_ptr<Tensor>> tensor721 = {I719, Gamma38_(), I720};
  auto task721 = make_shared<Task721>(tensor721, pindex);
  task720->add_dep(task721);
  task721->add_dep(task428);
  densityq->add_task(task721);

  vector<shared_ptr<Tensor>> tensor722 = {I720, t2};
  auto task722 = make_shared<Task722>(tensor722, pindex);
  task721->add_dep(task722);
  task722->add_dep(task428);
  densityq->add_task(task722);

  vector<IndexRange> I721_index = {virt_, virt_};
  auto I721 = make_shared<Tensor>(I721_index);
  vector<shared_ptr<Tensor>> tensor723 = {den2, I721};
  auto task723 = make_shared<Task723>(tensor723, pindex);
  task723->add_dep(task428);
  densityq->add_task(task723);

  vector<IndexRange> I722_index = {virt_, closed_, virt_, active_};
  auto I722 = make_shared<Tensor>(I722_index);
  vector<shared_ptr<Tensor>> tensor724 = {I721, t2, I722};
  auto task724 = make_shared<Task724>(tensor724, pindex);
  task723->add_dep(task724);
  task724->add_dep(task428);
  densityq->add_task(task724);

  vector<IndexRange> I723_index = {virt_, closed_, virt_, active_};
  auto I723 = make_shared<Tensor>(I723_index);
  vector<shared_ptr<Tensor>> tensor725 = {I722, Gamma38_(), I723};
  auto task725 = make_shared<Task725>(tensor725, pindex);
  task724->add_dep(task725);
  task725->add_dep(task428);
  densityq->add_task(task725);

  vector<shared_ptr<Tensor>> tensor726 = {I723, t2};
  auto task726 = make_shared<Task726>(tensor726, pindex);
  task725->add_dep(task726);
  task726->add_dep(task428);
  densityq->add_task(task726);

  vector<IndexRange> I725_index = {virt_, closed_, virt_, active_};
  auto I725 = make_shared<Tensor>(I725_index);
  vector<shared_ptr<Tensor>> tensor727 = {I721, t2, I725};
  auto task727 = make_shared<Task727>(tensor727, pindex);
  task723->add_dep(task727);
  task727->add_dep(task428);
  densityq->add_task(task727);

  vector<IndexRange> I726_index = {virt_, closed_, virt_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  vector<shared_ptr<Tensor>> tensor728 = {I725, Gamma38_(), I726};
  auto task728 = make_shared<Task728>(tensor728, pindex);
  task727->add_dep(task728);
  task728->add_dep(task428);
  densityq->add_task(task728);

  vector<shared_ptr<Tensor>> tensor729 = {I726, t2};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task728->add_dep(task729);
  task729->add_dep(task428);
  densityq->add_task(task729);

  vector<IndexRange> I727_index = {virt_, virt_};
  auto I727 = make_shared<Tensor>(I727_index);
  vector<shared_ptr<Tensor>> tensor730 = {den2, I727};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task730->add_dep(task428);
  densityq->add_task(task730);

  vector<IndexRange> I728_index = {virt_, closed_, virt_, active_};
  auto I728 = make_shared<Tensor>(I728_index);
  vector<shared_ptr<Tensor>> tensor731 = {I727, t2, I728};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task730->add_dep(task731);
  task731->add_dep(task428);
  densityq->add_task(task731);

  vector<IndexRange> I729_index = {virt_, closed_, virt_, active_};
  auto I729 = make_shared<Tensor>(I729_index);
  vector<shared_ptr<Tensor>> tensor732 = {I728, Gamma38_(), I729};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task731->add_dep(task732);
  task732->add_dep(task428);
  densityq->add_task(task732);

  vector<shared_ptr<Tensor>> tensor733 = {I729, t2};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  task733->add_dep(task428);
  densityq->add_task(task733);

  vector<IndexRange> I731_index = {virt_, closed_, virt_, active_};
  auto I731 = make_shared<Tensor>(I731_index);
  vector<shared_ptr<Tensor>> tensor734 = {I727, t2, I731};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task730->add_dep(task734);
  task734->add_dep(task428);
  densityq->add_task(task734);

  vector<IndexRange> I732_index = {virt_, closed_, virt_, active_};
  auto I732 = make_shared<Tensor>(I732_index);
  vector<shared_ptr<Tensor>> tensor735 = {I731, Gamma38_(), I732};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task734->add_dep(task735);
  task735->add_dep(task428);
  densityq->add_task(task735);

  vector<shared_ptr<Tensor>> tensor736 = {I732, t2};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task735->add_dep(task736);
  task736->add_dep(task428);
  densityq->add_task(task736);

  vector<IndexRange> I733_index = {active_, closed_};
  auto I733 = make_shared<Tensor>(I733_index);
  vector<shared_ptr<Tensor>> tensor737 = {den2, I733};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task737->add_dep(task428);
  densityq->add_task(task737);

  vector<IndexRange> I734_index = {virt_, virt_, active_, active_};
  auto I734 = make_shared<Tensor>(I734_index);
  vector<shared_ptr<Tensor>> tensor738 = {I733, t2, I734};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task737->add_dep(task738);
  task738->add_dep(task428);
  densityq->add_task(task738);

  vector<IndexRange> I735_index = {active_, virt_, active_, virt_};
  auto I735 = make_shared<Tensor>(I735_index);
  vector<shared_ptr<Tensor>> tensor739 = {I734, Gamma60_(), I735};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task738->add_dep(task739);
  task739->add_dep(task428);
  densityq->add_task(task739);

  vector<shared_ptr<Tensor>> tensor740 = {I735, t2};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task739->add_dep(task740);
  task740->add_dep(task428);
  densityq->add_task(task740);

  return densityq;
}


