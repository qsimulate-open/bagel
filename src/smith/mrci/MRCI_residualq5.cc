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
#include <src/smith/mrci/MRCI_tasks10.h>
#include <src/smith/mrci/MRCI_tasks11.h>
#include <src/smith/mrci/MRCI_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq5(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I108_index = {closed_, active_, active_, virt_};
  auto I108 = make_shared<Tensor>(I108_index);
  auto tensor451 = vector<shared_ptr<Tensor>>{r, I108};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task451->add_dep(task108);
  residualq->add_task(task451);

  vector<IndexRange> I109_index = {closed_, active_, active_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor452 = vector<shared_ptr<Tensor>>{I108, h1_, I109};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task451->add_dep(task452);
  task452->add_dep(task108);
  residualq->add_task(task452);

  auto tensor453 = vector<shared_ptr<Tensor>>{I109, Gamma4_(), t2};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task452->add_dep(task453);
  task453->add_dep(task108);
  residualq->add_task(task453);

  vector<IndexRange> I112_index = {active_, virt_, closed_, active_};
  auto I112 = make_shared<Tensor>(I112_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I108, Gamma5_(), I112};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task451->add_dep(task454);
  task454->add_dep(task108);
  residualq->add_task(task454);

  auto tensor455 = vector<shared_ptr<Tensor>>{I112, t2, h1_};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task108);
  residualq->add_task(task455);

  auto tensor456 = vector<shared_ptr<Tensor>>{I112, t2, h1_};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task454->add_dep(task456);
  task456->add_dep(task108);
  residualq->add_task(task456);

  auto tensor457 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task454->add_dep(task457);
  task457->add_dep(task108);
  residualq->add_task(task457);

  auto tensor458 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task454->add_dep(task458);
  task458->add_dep(task108);
  residualq->add_task(task458);

  auto tensor459 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task454->add_dep(task459);
  task459->add_dep(task108);
  residualq->add_task(task459);

  auto tensor460 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task454->add_dep(task460);
  task460->add_dep(task108);
  residualq->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task454->add_dep(task461);
  task461->add_dep(task108);
  residualq->add_task(task461);

  auto tensor462 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task454->add_dep(task462);
  task462->add_dep(task108);
  residualq->add_task(task462);

  auto tensor463 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task463 = make_shared<Task463>(tensor463, pindex);
  task454->add_dep(task463);
  task463->add_dep(task108);
  residualq->add_task(task463);

  auto tensor464 = vector<shared_ptr<Tensor>>{I112, t2, v2_};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task454->add_dep(task464);
  task464->add_dep(task108);
  residualq->add_task(task464);

  vector<IndexRange> I118_index = {closed_, active_, virt_, active_};
  auto I118 = make_shared<Tensor>(I118_index);
  auto tensor465 = vector<shared_ptr<Tensor>>{I108, Gamma29_(), I118};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task451->add_dep(task465);
  task465->add_dep(task108);
  residualq->add_task(task465);

  auto tensor466 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task465->add_dep(task466);
  task466->add_dep(task108);
  residualq->add_task(task466);

  auto tensor467 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task465->add_dep(task467);
  task467->add_dep(task108);
  residualq->add_task(task467);

  auto tensor468 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task465->add_dep(task468);
  task468->add_dep(task108);
  residualq->add_task(task468);

  auto tensor469 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task465->add_dep(task469);
  task469->add_dep(task108);
  residualq->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task465->add_dep(task470);
  task470->add_dep(task108);
  residualq->add_task(task470);

  auto tensor471 = vector<shared_ptr<Tensor>>{I118, t2, h1_};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task465->add_dep(task471);
  task471->add_dep(task108);
  residualq->add_task(task471);

  auto tensor472 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task465->add_dep(task472);
  task472->add_dep(task108);
  residualq->add_task(task472);

  auto tensor473 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task465->add_dep(task473);
  task473->add_dep(task108);
  residualq->add_task(task473);

  auto tensor474 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task465->add_dep(task474);
  task474->add_dep(task108);
  residualq->add_task(task474);

  auto tensor475 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task465->add_dep(task475);
  task475->add_dep(task108);
  residualq->add_task(task475);

  vector<IndexRange> I958_index = {virt_, closed_, active_, active_};
  auto I958 = make_shared<Tensor>(I958_index);
  auto tensor476 = vector<shared_ptr<Tensor>>{I118, t2, I958};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task465->add_dep(task476);
  task476->add_dep(task108);
  residualq->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I958, v2_};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task476->add_dep(task477);
  task477->add_dep(task108);
  residualq->add_task(task477);

  vector<IndexRange> I961_index = {virt_, closed_, active_, active_};
  auto I961 = make_shared<Tensor>(I961_index);
  auto tensor478 = vector<shared_ptr<Tensor>>{I118, t2, I961};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task465->add_dep(task478);
  task478->add_dep(task108);
  residualq->add_task(task478);

  auto tensor479 = vector<shared_ptr<Tensor>>{I961, v2_};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task478->add_dep(task479);
  task479->add_dep(task108);
  residualq->add_task(task479);

  auto tensor480 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task465->add_dep(task480);
  task480->add_dep(task108);
  residualq->add_task(task480);

  auto tensor481 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task465->add_dep(task481);
  task481->add_dep(task108);
  residualq->add_task(task481);

  vector<IndexRange> I1006_index = {virt_, active_, closed_, closed_};
  auto I1006 = make_shared<Tensor>(I1006_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{I118, t2, I1006};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task465->add_dep(task482);
  task482->add_dep(task108);
  residualq->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I1006, v2_};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task108);
  residualq->add_task(task483);

  vector<IndexRange> I1009_index = {virt_, active_, closed_, closed_};
  auto I1009 = make_shared<Tensor>(I1009_index);
  auto tensor484 = vector<shared_ptr<Tensor>>{I118, t2, I1009};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task465->add_dep(task484);
  task484->add_dep(task108);
  residualq->add_task(task484);

  auto tensor485 = vector<shared_ptr<Tensor>>{I1009, v2_};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task484->add_dep(task485);
  task485->add_dep(task108);
  residualq->add_task(task485);

  auto tensor486 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task465->add_dep(task486);
  task486->add_dep(task108);
  residualq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I118, t2, v2_};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task465->add_dep(task487);
  task487->add_dep(task108);
  residualq->add_task(task487);

  vector<IndexRange> I130_index = {virt_, active_, active_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor488 = vector<shared_ptr<Tensor>>{I108, h1_, I130};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task451->add_dep(task488);
  task488->add_dep(task108);
  residualq->add_task(task488);

  auto tensor489 = vector<shared_ptr<Tensor>>{I130, Gamma252_(), t2};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task488->add_dep(task489);
  task489->add_dep(task108);
  residualq->add_task(task489);

  vector<IndexRange> I133_index = {closed_, virt_};
  auto I133 = make_shared<Tensor>(I133_index);
  auto tensor490 = vector<shared_ptr<Tensor>>{I108, Gamma32_(), I133};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task451->add_dep(task490);
  task490->add_dep(task108);
  residualq->add_task(task490);

  auto tensor491 = vector<shared_ptr<Tensor>>{I133, t2, h1_};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task490->add_dep(task491);
  task491->add_dep(task108);
  residualq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I133, t2, h1_};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task490->add_dep(task492);
  task492->add_dep(task108);
  residualq->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I133, t2, v2_};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task490->add_dep(task493);
  task493->add_dep(task108);
  residualq->add_task(task493);

  auto tensor494 = vector<shared_ptr<Tensor>>{I133, t2, v2_};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task490->add_dep(task494);
  task494->add_dep(task108);
  residualq->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I133, t2, v2_};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task490->add_dep(task495);
  task495->add_dep(task108);
  residualq->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I133, t2, v2_};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task490->add_dep(task496);
  task496->add_dep(task108);
  residualq->add_task(task496);

  vector<IndexRange> I840_index = {closed_, closed_, active_, active_, active_, active_};
  auto I840 = make_shared<Tensor>(I840_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{I108, v2_, I840};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task451->add_dep(task497);
  task497->add_dep(task108);
  residualq->add_task(task497);

  auto tensor498 = vector<shared_ptr<Tensor>>{I840, Gamma107_(), t2};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task108);
  residualq->add_task(task498);

  vector<IndexRange> I843_index = {closed_, active_, active_, active_, active_, active_};
  auto I843 = make_shared<Tensor>(I843_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{I108, v2_, I843};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task451->add_dep(task499);
  task499->add_dep(task108);
  residualq->add_task(task499);

  auto tensor500 = vector<shared_ptr<Tensor>>{I843, Gamma278_(), t2};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task108);
  residualq->add_task(task500);

  vector<IndexRange> I846_index = {closed_, active_, active_, active_, active_, active_};
  auto I846 = make_shared<Tensor>(I846_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{I108, v2_, I846};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task451->add_dep(task501);
  task501->add_dep(task108);
  residualq->add_task(task501);

  auto tensor502 = vector<shared_ptr<Tensor>>{I846, Gamma100_(), t2};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task108);
  residualq->add_task(task502);

  vector<IndexRange> I849_index = {closed_, active_, active_, active_};
  auto I849 = make_shared<Tensor>(I849_index);
  auto tensor503 = vector<shared_ptr<Tensor>>{I108, v2_, I849};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task451->add_dep(task503);
  task503->add_dep(task108);
  residualq->add_task(task503);

  auto tensor504 = vector<shared_ptr<Tensor>>{I849, Gamma4_(), t2};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task503->add_dep(task504);
  task504->add_dep(task108);
  residualq->add_task(task504);

  vector<IndexRange> I852_index = {closed_, active_, active_, active_};
  auto I852 = make_shared<Tensor>(I852_index);
  auto tensor505 = vector<shared_ptr<Tensor>>{I108, v2_, I852};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task451->add_dep(task505);
  task505->add_dep(task108);
  residualq->add_task(task505);

  auto tensor506 = vector<shared_ptr<Tensor>>{I852, Gamma4_(), t2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task505->add_dep(task506);
  task506->add_dep(task108);
  residualq->add_task(task506);

  vector<IndexRange> I855_index = {closed_, active_, active_, active_};
  auto I855 = make_shared<Tensor>(I855_index);
  auto tensor507 = vector<shared_ptr<Tensor>>{I108, t2, I855};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task451->add_dep(task507);
  task507->add_dep(task108);
  residualq->add_task(task507);

  auto tensor508 = vector<shared_ptr<Tensor>>{I855, Gamma221_(), v2_};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task507->add_dep(task508);
  task508->add_dep(task108);
  residualq->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I855, Gamma104_(), v2_};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task507->add_dep(task509);
  task509->add_dep(task108);
  residualq->add_task(task509);

  vector<IndexRange> I858_index = {closed_, active_, active_, active_};
  auto I858 = make_shared<Tensor>(I858_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{I108, t2, I858};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task451->add_dep(task510);
  task510->add_dep(task108);
  residualq->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I858, Gamma221_(), v2_};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task510->add_dep(task511);
  task511->add_dep(task108);
  residualq->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I858, Gamma104_(), v2_};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task510->add_dep(task512);
  task512->add_dep(task108);
  residualq->add_task(task512);

  vector<IndexRange> I885_index = {closed_, closed_, active_, active_, active_, active_};
  auto I885 = make_shared<Tensor>(I885_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{I108, t2, I885};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task451->add_dep(task513);
  task513->add_dep(task108);
  residualq->add_task(task513);

  vector<IndexRange> I886_index = {closed_, closed_, active_, active_};
  auto I886 = make_shared<Tensor>(I886_index);
  auto tensor514 = vector<shared_ptr<Tensor>>{I885, Gamma240_(), I886};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task108);
  residualq->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I886, v2_};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task514->add_dep(task515);
  task515->add_dep(task108);
  residualq->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I885, Gamma7_(), v2_};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task513->add_dep(task516);
  task516->add_dep(task108);
  residualq->add_task(task516);

  auto tensor517 = vector<shared_ptr<Tensor>>{I885, Gamma296_(), v2_};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task513->add_dep(task517);
  task517->add_dep(task108);
  residualq->add_task(task517);

  vector<IndexRange> I888_index = {virt_, active_, active_, active_, closed_, active_};
  auto I888 = make_shared<Tensor>(I888_index);
  auto tensor518 = vector<shared_ptr<Tensor>>{I108, Gamma240_(), I888};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task451->add_dep(task518);
  task518->add_dep(task108);
  residualq->add_task(task518);

  vector<IndexRange> I889_index = {virt_, virt_, active_, active_};
  auto I889 = make_shared<Tensor>(I889_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{I888, t2, I889};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task518->add_dep(task519);
  task519->add_dep(task108);
  residualq->add_task(task519);

  auto tensor520 = vector<shared_ptr<Tensor>>{I889, v2_};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task519->add_dep(task520);
  task520->add_dep(task108);
  residualq->add_task(task520);

  vector<IndexRange> I919_index = {virt_, virt_, active_, active_};
  auto I919 = make_shared<Tensor>(I919_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{I888, t2, I919};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task518->add_dep(task521);
  task521->add_dep(task108);
  residualq->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I919, v2_};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task108);
  residualq->add_task(task522);

  vector<IndexRange> I894_index = {active_, active_, virt_, active_, closed_, active_};
  auto I894 = make_shared<Tensor>(I894_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{I108, Gamma7_(), I894};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task451->add_dep(task523);
  task523->add_dep(task108);
  residualq->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I894, t2, v2_};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task108);
  residualq->add_task(task524);

  vector<IndexRange> I900_index = {active_, virt_, active_, active_, closed_, active_};
  auto I900 = make_shared<Tensor>(I900_index);
  auto tensor525 = vector<shared_ptr<Tensor>>{I108, Gamma296_(), I900};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task451->add_dep(task525);
  task525->add_dep(task108);
  residualq->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I900, t2, v2_};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task525->add_dep(task526);
  task526->add_dep(task108);
  residualq->add_task(task526);

  vector<IndexRange> I915_index = {closed_, closed_, active_, active_, active_, active_};
  auto I915 = make_shared<Tensor>(I915_index);
  auto tensor527 = vector<shared_ptr<Tensor>>{I108, t2, I915};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task451->add_dep(task527);
  task527->add_dep(task108);
  residualq->add_task(task527);

  vector<IndexRange> I916_index = {closed_, closed_, active_, active_};
  auto I916 = make_shared<Tensor>(I916_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{I915, Gamma240_(), I916};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task108);
  residualq->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I916, v2_};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task528->add_dep(task529);
  task529->add_dep(task108);
  residualq->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I915, Gamma4_(), v2_};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task527->add_dep(task530);
  task530->add_dep(task108);
  residualq->add_task(task530);

  vector<IndexRange> I924_index = {active_, active_, virt_, closed_, active_, active_};
  auto I924 = make_shared<Tensor>(I924_index);
  auto tensor531 = vector<shared_ptr<Tensor>>{I108, Gamma4_(), I924};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task451->add_dep(task531);
  task531->add_dep(task108);
  residualq->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I924, t2, v2_};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task531->add_dep(task532);
  task532->add_dep(task108);
  residualq->add_task(task532);

  vector<IndexRange> I945_index = {closed_, active_, active_, active_, active_, active_};
  auto I945 = make_shared<Tensor>(I945_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{I108, t2, I945};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task451->add_dep(task533);
  task533->add_dep(task108);
  residualq->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I945, Gamma312_(), v2_};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task108);
  residualq->add_task(task534);

  auto tensor535 = vector<shared_ptr<Tensor>>{I945, Gamma313_(), v2_};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task533->add_dep(task535);
  task535->add_dep(task108);
  residualq->add_task(task535);

  vector<IndexRange> I951_index = {virt_, active_, active_, active_};
  auto I951 = make_shared<Tensor>(I951_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{I108, v2_, I951};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task451->add_dep(task536);
  task536->add_dep(task108);
  residualq->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I951, Gamma252_(), t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task108);
  residualq->add_task(task537);

  vector<IndexRange> I954_index = {virt_, active_, active_, active_};
  auto I954 = make_shared<Tensor>(I954_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{I108, v2_, I954};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task451->add_dep(task538);
  task538->add_dep(task108);
  residualq->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I954, Gamma252_(), t2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task108);
  residualq->add_task(task539);

  vector<IndexRange> I993_index = {virt_, active_, active_, active_};
  auto I993 = make_shared<Tensor>(I993_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{I108, t2, I993};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task451->add_dep(task540);
  task540->add_dep(task108);
  residualq->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I993, Gamma240_(), v2_};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task108);
  residualq->add_task(task541);

  auto tensor542 = vector<shared_ptr<Tensor>>{I993, Gamma252_(), v2_};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task540->add_dep(task542);
  task542->add_dep(task108);
  residualq->add_task(task542);

  vector<IndexRange> I996_index = {virt_, active_, active_, active_};
  auto I996 = make_shared<Tensor>(I996_index);
  auto tensor543 = vector<shared_ptr<Tensor>>{I108, t2, I996};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task451->add_dep(task543);
  task543->add_dep(task108);
  residualq->add_task(task543);

  auto tensor544 = vector<shared_ptr<Tensor>>{I996, Gamma240_(), v2_};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task543->add_dep(task544);
  task544->add_dep(task108);
  residualq->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I996, Gamma252_(), v2_};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task543->add_dep(task545);
  task545->add_dep(task108);
  residualq->add_task(task545);

  vector<IndexRange> I1023_index = {closed_, active_, active_, active_, virt_, active_};
  auto I1023 = make_shared<Tensor>(I1023_index);
  auto tensor546 = vector<shared_ptr<Tensor>>{I108, Gamma338_(), I1023};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task451->add_dep(task546);
  task546->add_dep(task108);
  residualq->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I1023, t2, v2_};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task546->add_dep(task547);
  task547->add_dep(task108);
  residualq->add_task(task547);

  vector<IndexRange> I1717_index = {active_, virt_, closed_, active_};
  auto I1717 = make_shared<Tensor>(I1717_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{I108, Gamma569_(), I1717};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task451->add_dep(task548);
  task548->add_dep(task108);
  residualq->add_task(task548);

  auto tensor549 = vector<shared_ptr<Tensor>>{I1717, t2};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task548->add_dep(task549);
  task549->add_dep(task108);
  residualq->add_task(task549);

  vector<IndexRange> I1725_index = {active_, virt_, closed_, active_};
  auto I1725 = make_shared<Tensor>(I1725_index);
  auto tensor550 = vector<shared_ptr<Tensor>>{I108, Gamma573_(), I1725};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task451->add_dep(task550);
  task550->add_dep(task108);
  residualq->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I1725, t2};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task550->add_dep(task551);
  task551->add_dep(task108);
  residualq->add_task(task551);
}

#endif
