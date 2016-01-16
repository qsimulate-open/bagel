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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_deciq(const bool reset, const bool diagonal) {

  auto deciq = make_shared<Queue>();
  auto task537 = make_shared<Task537>(deci, reset);
  deciq->add_task(task537);

  auto I698 = make_shared<TATensor<double,1>>({ci_});
  auto task538 = make_shared<Task538>(deci, I698);
  task538->add_dep(task537);
  deciq->add_task(task538);

  auto I699 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task539 = make_shared<Task539>(I698, Gamma248_(), I699);
  task538->add_dep(task539);
  task539->add_dep(task537);
  deciq->add_task(task539);

  auto task540 = make_shared<Task540>(I699, t2);
  task539->add_dep(task540);
  task540->add_dep(task537);
  deciq->add_task(task540);

  auto I702 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task541 = make_shared<Task541>(I698, Gamma249_(), I702);
  task538->add_dep(task541);
  task541->add_dep(task537);
  deciq->add_task(task541);

  auto I703 = make_shared<TATensor<double,4>>({active_, active_, closed_, closed_});
  auto task542 = make_shared<Task542>(I702, t2, I703);
  task541->add_dep(task542);
  task542->add_dep(task537);
  deciq->add_task(task542);

  auto task543 = make_shared<Task543>(I703, f1_, t2);
  task542->add_dep(task543);
  task543->add_dep(task537);
  deciq->add_task(task543);

  auto task544 = make_shared<Task544>(I702, t2, this->e0_);
  task541->add_dep(task544);
  task544->add_dep(task537);
  deciq->add_task(task544);

  auto task545 = make_shared<Task545>(I702, v2_, t2);
  task541->add_dep(task545);
  task545->add_dep(task537);
  deciq->add_task(task545);

  auto task546 = make_shared<Task546>(I702, v2_, t2);
  task541->add_dep(task546);
  task546->add_dep(task537);
  deciq->add_task(task546);

  auto I706 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task547 = make_shared<Task547>(I698, Gamma250_(), I706);
  task538->add_dep(task547);
  task547->add_dep(task537);
  deciq->add_task(task547);

  auto I707 = make_shared<TATensor<double,4>>({active_, active_, closed_, active_});
  auto task548 = make_shared<Task548>(I706, t2, I707);
  task547->add_dep(task548);
  task548->add_dep(task537);
  deciq->add_task(task548);

  auto task549 = make_shared<Task549>(I707, f1_, t2);
  task548->add_dep(task549);
  task549->add_dep(task537);
  deciq->add_task(task549);

  auto I710 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task550 = make_shared<Task550>(I698, Gamma251_(), I710);
  task538->add_dep(task550);
  task550->add_dep(task537);
  deciq->add_task(task550);

  auto I711 = make_shared<TATensor<double,4>>({active_, closed_, closed_, active_});
  auto task551 = make_shared<Task551>(I710, t2, I711);
  task550->add_dep(task551);
  task551->add_dep(task537);
  deciq->add_task(task551);

  auto task552 = make_shared<Task552>(I711, t2, f1_);
  task551->add_dep(task552);
  task552->add_dep(task537);
  deciq->add_task(task552);

  auto I742 = make_shared<TATensor<double,4>>({active_, closed_, closed_, active_});
  auto task553 = make_shared<Task553>(I710, t2, I742);
  task550->add_dep(task553);
  task553->add_dep(task537);
  deciq->add_task(task553);

  auto task554 = make_shared<Task554>(I742, f1_, t2);
  task553->add_dep(task554);
  task554->add_dep(task537);
  deciq->add_task(task554);

  auto I714 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task555 = make_shared<Task555>(I698, Gamma252_(), I714);
  task538->add_dep(task555);
  task555->add_dep(task537);
  deciq->add_task(task555);

  auto I715 = make_shared<TATensor<double,4>>({active_, closed_, active_, active_});
  auto task556 = make_shared<Task556>(I714, t2, I715);
  task555->add_dep(task556);
  task556->add_dep(task537);
  deciq->add_task(task556);

  auto task557 = make_shared<Task557>(I715, t2, f1_);
  task556->add_dep(task557);
  task557->add_dep(task537);
  deciq->add_task(task557);

  auto I718 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task558 = make_shared<Task558>(I698, Gamma253_(), I718);
  task538->add_dep(task558);
  task558->add_dep(task537);
  deciq->add_task(task558);

  auto task559 = make_shared<Task559>(I718, t2);
  task558->add_dep(task559);
  task559->add_dep(task537);
  deciq->add_task(task559);

  auto I721 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task560 = make_shared<Task560>(I698, Gamma254_(), I721);
  task538->add_dep(task560);
  task560->add_dep(task537);
  deciq->add_task(task560);

  auto I722 = make_shared<TATensor<double,4>>({active_, active_, active_, closed_});
  auto task561 = make_shared<Task561>(I721, t2, I722);
  task560->add_dep(task561);
  task561->add_dep(task537);
  deciq->add_task(task561);

  auto task562 = make_shared<Task562>(I722, f1_, t2);
  task561->add_dep(task562);
  task562->add_dep(task537);
  deciq->add_task(task562);

  auto I738 = make_shared<TATensor<double,4>>({active_, closed_, active_, active_});
  auto task563 = make_shared<Task563>(I721, t2, I738);
  task560->add_dep(task563);
  task563->add_dep(task537);
  deciq->add_task(task563);

  auto task564 = make_shared<Task564>(I738, t2, f1_);
  task563->add_dep(task564);
  task564->add_dep(task537);
  deciq->add_task(task564);

  auto I862 = make_shared<TATensor<double,4>>({active_, active_, closed_, active_});
  auto task565 = make_shared<Task565>(I721, t2, I862);
  task560->add_dep(task565);
  task565->add_dep(task537);
  deciq->add_task(task565);

  auto task566 = make_shared<Task566>(I862, t2, this->e0_);
  task565->add_dep(task566);
  task566->add_dep(task537);
  deciq->add_task(task566);

  auto task567 = make_shared<Task567>(I862, f1_, t2);
  task565->add_dep(task567);
  task567->add_dep(task537);
  deciq->add_task(task567);

  auto task568 = make_shared<Task568>(I721, v2_, t2);
  task560->add_dep(task568);
  task568->add_dep(task537);
  deciq->add_task(task568);

  auto task569 = make_shared<Task569>(I721, v2_, t2);
  task560->add_dep(task569);
  task569->add_dep(task537);
  deciq->add_task(task569);

  auto I725 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task570 = make_shared<Task570>(I698, Gamma255_(), I725);
  task538->add_dep(task570);
  task570->add_dep(task537);
  deciq->add_task(task570);

  auto I726 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task571 = make_shared<Task571>(I725, t2, I726);
  task570->add_dep(task571);
  task571->add_dep(task537);
  deciq->add_task(task571);

  auto task572 = make_shared<Task572>(I726, t2, f1_);
  task571->add_dep(task572);
  task572->add_dep(task537);
  deciq->add_task(task572);

  auto task573 = make_shared<Task573>(I726, t2, f1_);
  task571->add_dep(task573);
  task573->add_dep(task537);
  deciq->add_task(task573);

  auto I816 = make_shared<TATensor<double,4>>({active_, closed_, virt_, active_});
  auto task574 = make_shared<Task574>(I725, t2, I816);
  task570->add_dep(task574);
  task574->add_dep(task537);
  deciq->add_task(task574);

  auto task575 = make_shared<Task575>(I816, t2, f1_);
  task574->add_dep(task575);
  task575->add_dep(task537);
  deciq->add_task(task575);

  auto I866 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task576 = make_shared<Task576>(I725, t2, I866);
  task570->add_dep(task576);
  task576->add_dep(task537);
  deciq->add_task(task576);

  auto task577 = make_shared<Task577>(I866, t2, f1_);
  task576->add_dep(task577);
  task577->add_dep(task537);
  deciq->add_task(task577);

  auto task578 = make_shared<Task578>(I866, t2, f1_);
  task576->add_dep(task578);
  task578->add_dep(task537);
  deciq->add_task(task578);

  auto task579 = make_shared<Task579>(I725, v2_, t2);
  task570->add_dep(task579);
  task579->add_dep(task537);
  deciq->add_task(task579);

  auto task580 = make_shared<Task580>(I725, h1_, t2);
  task570->add_dep(task580);
  task580->add_dep(task537);
  deciq->add_task(task580);

  auto I733 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task581 = make_shared<Task581>(I698, Gamma257_(), I733);
  task538->add_dep(task581);
  task581->add_dep(task537);
  deciq->add_task(task581);

  auto I734 = make_shared<TATensor<double,4>>({active_, active_, closed_, active_});
  auto task582 = make_shared<Task582>(I733, t2, I734);
  task581->add_dep(task582);
  task582->add_dep(task537);
  deciq->add_task(task582);

  auto task583 = make_shared<Task583>(I734, t2, f1_);
  task582->add_dep(task583);
  task583->add_dep(task537);
  deciq->add_task(task583);

  auto I745 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task584 = make_shared<Task584>(I698, Gamma260_(), I745);
  task538->add_dep(task584);
  task584->add_dep(task537);
  deciq->add_task(task584);

  auto I746 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task585 = make_shared<Task585>(I745, t2, I746);
  task584->add_dep(task585);
  task585->add_dep(task537);
  deciq->add_task(task585);

  auto task586 = make_shared<Task586>(I746, f1_, t2);
  task585->add_dep(task586);
  task586->add_dep(task537);
  deciq->add_task(task586);

  auto I750 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task587 = make_shared<Task587>(I745, t2, I750);
  task584->add_dep(task587);
  task587->add_dep(task537);
  deciq->add_task(task587);

  auto task588 = make_shared<Task588>(I750, f1_, t2);
  task587->add_dep(task588);
  task588->add_dep(task537);
  deciq->add_task(task588);

  auto I788 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task589 = make_shared<Task589>(I745, t2, I788);
  task584->add_dep(task589);
  task589->add_dep(task537);
  deciq->add_task(task589);

  auto task590 = make_shared<Task590>(I788, f1_, t2);
  task589->add_dep(task590);
  task590->add_dep(task537);
  deciq->add_task(task590);

  auto I792 = make_shared<TATensor<double,4>>({active_, closed_, virt_, active_});
  auto task591 = make_shared<Task591>(I745, t2, I792);
  task584->add_dep(task591);
  task591->add_dep(task537);
  deciq->add_task(task591);

  auto task592 = make_shared<Task592>(I792, f1_, t2);
  task591->add_dep(task592);
  task592->add_dep(task537);
  deciq->add_task(task592);

  auto I796 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task593 = make_shared<Task593>(I745, t2, I796);
  task584->add_dep(task593);
  task593->add_dep(task537);
  deciq->add_task(task593);

  auto task594 = make_shared<Task594>(I796, f1_, t2);
  task593->add_dep(task594);
  task594->add_dep(task537);
  deciq->add_task(task594);

  auto task595 = make_shared<Task595>(I745, v2_, t2);
  task584->add_dep(task595);
  task595->add_dep(task537);
  deciq->add_task(task595);

  auto task596 = make_shared<Task596>(I745, h1_, t2);
  task584->add_dep(task596);
  task596->add_dep(task537);
  deciq->add_task(task596);

  auto I753 = make_shared<TATensor<double,2>>({active_, active_});
  auto task597 = make_shared<Task597>(I698, Gamma262_(), I753);
  task538->add_dep(task597);
  task597->add_dep(task537);
  deciq->add_task(task597);

  auto task598 = make_shared<Task598>(I753, t2);
  task597->add_dep(task598);
  task598->add_dep(task537);
  deciq->add_task(task598);

  auto task599 = make_shared<Task599>(I753, t2);
  task597->add_dep(task599);
  task599->add_dep(task537);
  deciq->add_task(task599);

  auto I759 = make_shared<TATensor<double,2>>({active_, active_});
  auto task600 = make_shared<Task600>(I698, Gamma264_(), I759);
  task538->add_dep(task600);
  task600->add_dep(task537);
  deciq->add_task(task600);

  auto I760 = make_shared<TATensor<double,4>>({active_, closed_, virt_, closed_});
  auto task601 = make_shared<Task601>(I759, t2, I760);
  task600->add_dep(task601);
  task601->add_dep(task537);
  deciq->add_task(task601);

  auto task602 = make_shared<Task602>(I760, f1_, t2);
  task601->add_dep(task602);
  task602->add_dep(task537);
  deciq->add_task(task602);

  auto I764 = make_shared<TATensor<double,4>>({active_, closed_, virt_, closed_});
  auto task603 = make_shared<Task603>(I759, t2, I764);
  task600->add_dep(task603);
  task603->add_dep(task537);
  deciq->add_task(task603);

  auto task604 = make_shared<Task604>(I764, f1_, t2);
  task603->add_dep(task604);
  task604->add_dep(task537);
  deciq->add_task(task604);

  auto I768 = make_shared<TATensor<double,4>>({active_, virt_, closed_, closed_});
  auto task605 = make_shared<Task605>(I759, t2, I768);
  task600->add_dep(task605);
  task605->add_dep(task537);
  deciq->add_task(task605);

  auto task606 = make_shared<Task606>(I768, f1_, t2);
  task605->add_dep(task606);
  task606->add_dep(task537);
  deciq->add_task(task606);

  auto I772 = make_shared<TATensor<double,4>>({active_, closed_, closed_, virt_});
  auto task607 = make_shared<Task607>(I759, t2, I772);
  task600->add_dep(task607);
  task607->add_dep(task537);
  deciq->add_task(task607);

  auto task608 = make_shared<Task608>(I772, f1_, t2);
  task607->add_dep(task608);
  task608->add_dep(task537);
  deciq->add_task(task608);

  auto I776 = make_shared<TATensor<double,4>>({active_, virt_, closed_, closed_});
  auto task609 = make_shared<Task609>(I759, t2, I776);
  task600->add_dep(task609);
  task609->add_dep(task537);
  deciq->add_task(task609);

  auto task610 = make_shared<Task610>(I776, f1_, t2);
  task609->add_dep(task610);
  task610->add_dep(task537);
  deciq->add_task(task610);

  auto I780 = make_shared<TATensor<double,4>>({active_, closed_, closed_, virt_});
  auto task611 = make_shared<Task611>(I759, t2, I780);
  task600->add_dep(task611);
  task611->add_dep(task537);
  deciq->add_task(task611);

  auto task612 = make_shared<Task612>(I780, f1_, t2);
  task611->add_dep(task612);
  task612->add_dep(task537);
  deciq->add_task(task612);

  auto I800 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task613 = make_shared<Task613>(I759, f1_, I800);
  task600->add_dep(task613);
  task613->add_dep(task537);
  deciq->add_task(task613);

  auto task614 = make_shared<Task614>(I800, t2);
  task613->add_dep(task614);
  task614->add_dep(task537);
  deciq->add_task(task614);

  auto task615 = make_shared<Task615>(I800, t2);
  task613->add_dep(task615);
  task615->add_dep(task537);
  deciq->add_task(task615);

  auto I943 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task616 = make_shared<Task616>(I759, f1_, I943);
  task600->add_dep(task616);
  task616->add_dep(task537);
  deciq->add_task(task616);

  auto task617 = make_shared<Task617>(I943, t2);
  task616->add_dep(task617);
  task617->add_dep(task537);
  deciq->add_task(task617);

  auto I947 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task618 = make_shared<Task618>(I759, f1_, I947);
  task600->add_dep(task618);
  task618->add_dep(task537);
  deciq->add_task(task618);

  auto task619 = make_shared<Task619>(I947, t2);
  task618->add_dep(task619);
  task619->add_dep(task537);
  deciq->add_task(task619);

  auto task620 = make_shared<Task620>(I759, t2, this->e0_);
  task600->add_dep(task620);
  task620->add_dep(task537);
  deciq->add_task(task620);

  auto task621 = make_shared<Task621>(I759, t2, this->e0_);
  task600->add_dep(task621);
  task621->add_dep(task537);
  deciq->add_task(task621);

  auto task622 = make_shared<Task622>(I759, v2_, t2);
  task600->add_dep(task622);
  task622->add_dep(task537);
  deciq->add_task(task622);

  auto task623 = make_shared<Task623>(I759, v2_, t2);
  task600->add_dep(task623);
  task623->add_dep(task537);
  deciq->add_task(task623);

  auto task624 = make_shared<Task624>(I759, v2_, t2);
  task600->add_dep(task624);
  task624->add_dep(task537);
  deciq->add_task(task624);

  auto task625 = make_shared<Task625>(I759, v2_, t2);
  task600->add_dep(task625);
  task625->add_dep(task537);
  deciq->add_task(task625);

  auto I783 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task626 = make_shared<Task626>(I698, Gamma270_(), I783);
  task538->add_dep(task626);
  task626->add_dep(task537);
  deciq->add_task(task626);

  auto I784 = make_shared<TATensor<double,4>>({active_, closed_, virt_, active_});
  auto task627 = make_shared<Task627>(I783, t2, I784);
  task626->add_dep(task627);
  task627->add_dep(task537);
  deciq->add_task(task627);

  auto task628 = make_shared<Task628>(I784, f1_, t2);
  task627->add_dep(task628);
  task628->add_dep(task537);
  deciq->add_task(task628);

  auto task629 = make_shared<Task629>(I783, v2_, t2);
  task626->add_dep(task629);
  task629->add_dep(task537);
  deciq->add_task(task629);

  auto I807 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task630 = make_shared<Task630>(I698, Gamma276_(), I807);
  task538->add_dep(task630);
  task630->add_dep(task537);
  deciq->add_task(task630);

  auto I808 = make_shared<TATensor<double,4>>({active_, closed_, active_, active_});
  auto task631 = make_shared<Task631>(I807, t2, I808);
  task630->add_dep(task631);
  task631->add_dep(task537);
  deciq->add_task(task631);

  auto task632 = make_shared<Task632>(I808, f1_, t2);
  task631->add_dep(task632);
  task632->add_dep(task537);
  deciq->add_task(task632);

  auto I811 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task633 = make_shared<Task633>(I698, Gamma277_(), I811);
  task538->add_dep(task633);
  task633->add_dep(task537);
  deciq->add_task(task633);

  auto I812 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task634 = make_shared<Task634>(I811, t2, I812);
  task633->add_dep(task634);
  task634->add_dep(task537);
  deciq->add_task(task634);

  auto task635 = make_shared<Task635>(I812, t2, f1_);
  task634->add_dep(task635);
  task635->add_dep(task537);
  deciq->add_task(task635);

  auto task636 = make_shared<Task636>(I811, v2_, t2);
  task633->add_dep(task636);
  task636->add_dep(task537);
  deciq->add_task(task636);

  auto I819 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task637 = make_shared<Task637>(I698, Gamma279_(), I819);
  task538->add_dep(task637);
  task637->add_dep(task537);
  deciq->add_task(task637);

  auto task638 = make_shared<Task638>(I819, t2);
  task637->add_dep(task638);
  task638->add_dep(task537);
  deciq->add_task(task638);

  auto I822 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task639 = make_shared<Task639>(I698, Gamma280_(), I822);
  task538->add_dep(task639);
  task639->add_dep(task537);
  deciq->add_task(task639);

  auto I823 = make_shared<TATensor<double,4>>({active_, virt_, active_, closed_});
  auto task640 = make_shared<Task640>(I822, t2, I823);
  task639->add_dep(task640);
  task640->add_dep(task537);
  deciq->add_task(task640);

  auto task641 = make_shared<Task641>(I823, f1_, t2);
  task640->add_dep(task641);
  task641->add_dep(task537);
  deciq->add_task(task641);

  auto I827 = make_shared<TATensor<double,4>>({active_, closed_, active_, virt_});
  auto task642 = make_shared<Task642>(I822, t2, I827);
  task639->add_dep(task642);
  task642->add_dep(task537);
  deciq->add_task(task642);

  auto task643 = make_shared<Task643>(I827, f1_, t2);
  task642->add_dep(task643);
  task643->add_dep(task537);
  deciq->add_task(task643);

  auto I858 = make_shared<TATensor<double,4>>({active_, active_, virt_, closed_});
  auto task644 = make_shared<Task644>(I822, t2, I858);
  task639->add_dep(task644);
  task644->add_dep(task537);
  deciq->add_task(task644);

  auto task645 = make_shared<Task645>(I858, t2, f1_);
  task644->add_dep(task645);
  task645->add_dep(task537);
  deciq->add_task(task645);

  auto I985 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task646 = make_shared<Task646>(I822, t2, I985);
  task639->add_dep(task646);
  task646->add_dep(task537);
  deciq->add_task(task646);

  auto task647 = make_shared<Task647>(I985, t2, this->e0_);
  task646->add_dep(task647);
  task647->add_dep(task537);
  deciq->add_task(task647);

  auto task648 = make_shared<Task648>(I985, f1_, t2);
  task646->add_dep(task648);
  task648->add_dep(task537);
  deciq->add_task(task648);

  auto task649 = make_shared<Task649>(I822, v2_, t2);
  task639->add_dep(task649);
  task649->add_dep(task537);
  deciq->add_task(task649);

  auto task650 = make_shared<Task650>(I822, v2_, t2);
  task639->add_dep(task650);
  task650->add_dep(task537);
  deciq->add_task(task650);

  auto I830 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task651 = make_shared<Task651>(I698, Gamma282_(), I830);
  task538->add_dep(task651);
  task651->add_dep(task537);
  deciq->add_task(task651);

  auto task652 = make_shared<Task652>(I830, t2);
  task651->add_dep(task652);
  task652->add_dep(task537);
  deciq->add_task(task652);

  auto task653 = make_shared<Task653>(I830, t2);
  task651->add_dep(task653);
  task653->add_dep(task537);
  deciq->add_task(task653);

  auto task654 = make_shared<Task654>(I830, t2);
  task651->add_dep(task654);
  task654->add_dep(task537);
  deciq->add_task(task654);

  auto I833 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task655 = make_shared<Task655>(I698, Gamma283_(), I833);
  task538->add_dep(task655);
  task655->add_dep(task537);
  deciq->add_task(task655);

  auto I834 = make_shared<TATensor<double,4>>({active_, virt_, active_, closed_});
  auto task656 = make_shared<Task656>(I833, t2, I834);
  task655->add_dep(task656);
  task656->add_dep(task537);
  deciq->add_task(task656);

  auto task657 = make_shared<Task657>(I834, f1_, t2);
  task656->add_dep(task657);
  task657->add_dep(task537);
  deciq->add_task(task657);

  auto I838 = make_shared<TATensor<double,4>>({active_, closed_, active_, virt_});
  auto task658 = make_shared<Task658>(I833, t2, I838);
  task655->add_dep(task658);
  task658->add_dep(task537);
  deciq->add_task(task658);

  auto task659 = make_shared<Task659>(I838, f1_, t2);
  task658->add_dep(task659);
  task659->add_dep(task537);
  deciq->add_task(task659);

  auto task660 = make_shared<Task660>(I838, f1_, t2);
  task658->add_dep(task660);
  task660->add_dep(task537);
  deciq->add_task(task660);

  auto I854 = make_shared<TATensor<double,4>>({active_, active_, closed_, virt_});
  auto task661 = make_shared<Task661>(I833, t2, I854);
  task655->add_dep(task661);
  task661->add_dep(task537);
  deciq->add_task(task661);

  auto task662 = make_shared<Task662>(I854, t2, f1_);
  task661->add_dep(task662);
  task662->add_dep(task537);
  deciq->add_task(task662);

  auto I877 = make_shared<TATensor<double,4>>({active_, active_, virt_, closed_});
  auto task663 = make_shared<Task663>(I833, t2, I877);
  task655->add_dep(task663);
  task663->add_dep(task537);
  deciq->add_task(task663);

  auto task664 = make_shared<Task664>(I877, f1_, t2);
  task663->add_dep(task664);
  task664->add_dep(task537);
  deciq->add_task(task664);

  auto I881 = make_shared<TATensor<double,4>>({active_, active_, closed_, virt_});
  auto task665 = make_shared<Task665>(I833, t2, I881);
  task655->add_dep(task665);
  task665->add_dep(task537);
  deciq->add_task(task665);

  auto task666 = make_shared<Task666>(I881, f1_, t2);
  task665->add_dep(task666);
  task666->add_dep(task537);
  deciq->add_task(task666);

  auto I888 = make_shared<TATensor<double,4>>({active_, active_, virt_, closed_});
  auto task667 = make_shared<Task667>(I833, t2, I888);
  task655->add_dep(task667);
  task667->add_dep(task537);
  deciq->add_task(task667);

  auto task668 = make_shared<Task668>(I888, f1_, t2);
  task667->add_dep(task668);
  task668->add_dep(task537);
  deciq->add_task(task668);

  auto I892 = make_shared<TATensor<double,4>>({active_, active_, closed_, virt_});
  auto task669 = make_shared<Task669>(I833, t2, I892);
  task655->add_dep(task669);
  task669->add_dep(task537);
  deciq->add_task(task669);

  auto task670 = make_shared<Task670>(I892, f1_, t2);
  task669->add_dep(task670);
  task670->add_dep(task537);
  deciq->add_task(task670);

  auto I908 = make_shared<TATensor<double,4>>({active_, active_, closed_, virt_});
  auto task671 = make_shared<Task671>(I833, t2, I908);
  task655->add_dep(task671);
  task671->add_dep(task537);
  deciq->add_task(task671);

  auto task672 = make_shared<Task672>(I908, t2, f1_);
  task671->add_dep(task672);
  task672->add_dep(task537);
  deciq->add_task(task672);

  auto task673 = make_shared<Task673>(I908, t2, f1_);
  task671->add_dep(task673);
  task673->add_dep(task537);
  deciq->add_task(task673);

  auto I981 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task674 = make_shared<Task674>(I833, t2, I981);
  task655->add_dep(task674);
  task674->add_dep(task537);
  deciq->add_task(task674);

  auto task675 = make_shared<Task675>(I981, f1_, t2);
  task674->add_dep(task675);
  task675->add_dep(task537);
  deciq->add_task(task675);

  auto I993 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task676 = make_shared<Task676>(I833, t2, I993);
  task655->add_dep(task676);
  task676->add_dep(task537);
  deciq->add_task(task676);

  auto task677 = make_shared<Task677>(I993, t2, this->e0_);
  task676->add_dep(task677);
  task677->add_dep(task537);
  deciq->add_task(task677);

  auto task678 = make_shared<Task678>(I993, f1_, t2);
  task676->add_dep(task678);
  task678->add_dep(task537);
  deciq->add_task(task678);

  auto task679 = make_shared<Task679>(I833, t2, this->e0_);
  task655->add_dep(task679);
  task679->add_dep(task537);
  deciq->add_task(task679);

  auto task680 = make_shared<Task680>(I833, t2, this->e0_);
  task655->add_dep(task680);
  task680->add_dep(task537);
  deciq->add_task(task680);

  auto task681 = make_shared<Task681>(I833, v2_, t2);
  task655->add_dep(task681);
  task681->add_dep(task537);
  deciq->add_task(task681);

  auto task682 = make_shared<Task682>(I833, v2_, t2);
  task655->add_dep(task682);
  task682->add_dep(task537);
  deciq->add_task(task682);

  auto task683 = make_shared<Task683>(I833, v2_, t2);
  task655->add_dep(task683);
  task683->add_dep(task537);
  deciq->add_task(task683);

  auto task684 = make_shared<Task684>(I833, v2_, t2);
  task655->add_dep(task684);
  task684->add_dep(task537);
  deciq->add_task(task684);

  auto task685 = make_shared<Task685>(I833, v2_, t2);
  task655->add_dep(task685);
  task685->add_dep(task537);
  deciq->add_task(task685);

  auto task686 = make_shared<Task686>(I833, v2_, t2);
  task655->add_dep(task686);
  task686->add_dep(task537);
  deciq->add_task(task686);

  auto task687 = make_shared<Task687>(I833, v2_, t2);
  task655->add_dep(task687);
  task687->add_dep(task537);
  deciq->add_task(task687);

  auto task688 = make_shared<Task688>(I833, v2_, t2);
  task655->add_dep(task688);
  task688->add_dep(task537);
  deciq->add_task(task688);

  auto task689 = make_shared<Task689>(I833, v2_, t2);
  task655->add_dep(task689);
  task689->add_dep(task537);
  deciq->add_task(task689);

  auto task690 = make_shared<Task690>(I833, v2_, t2);
  task655->add_dep(task690);
  task690->add_dep(task537);
  deciq->add_task(task690);

  auto I841 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task691 = make_shared<Task691>(I698, Gamma285_(), I841);
  task538->add_dep(task691);
  task691->add_dep(task537);
  deciq->add_task(task691);

  auto I842 = make_shared<TATensor<double,4>>({active_, virt_, active_, active_});
  auto task692 = make_shared<Task692>(I841, t2, I842);
  task691->add_dep(task692);
  task692->add_dep(task537);
  deciq->add_task(task692);

  auto task693 = make_shared<Task693>(I842, f1_, t2);
  task692->add_dep(task693);
  task693->add_dep(task537);
  deciq->add_task(task693);

  auto I845 = make_shared<TATensor<double,2>>({active_, active_});
  auto task694 = make_shared<Task694>(I698, Gamma286_(), I845);
  task538->add_dep(task694);
  task694->add_dep(task537);
  deciq->add_task(task694);

  auto I846 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task695 = make_shared<Task695>(I845, t2, I846);
  task694->add_dep(task695);
  task695->add_dep(task537);
  deciq->add_task(task695);

  auto task696 = make_shared<Task696>(I846, t2, f1_);
  task695->add_dep(task696);
  task696->add_dep(task537);
  deciq->add_task(task696);

  auto task697 = make_shared<Task697>(I846, t2, f1_);
  task695->add_dep(task697);
  task697->add_dep(task537);
  deciq->add_task(task697);

  auto I900 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task698 = make_shared<Task698>(I845, t2, I900);
  task694->add_dep(task698);
  task698->add_dep(task537);
  deciq->add_task(task698);

  auto task699 = make_shared<Task699>(I900, t2, f1_);
  task698->add_dep(task699);
  task699->add_dep(task537);
  deciq->add_task(task699);

  auto task700 = make_shared<Task700>(I900, t2, f1_);
  task698->add_dep(task700);
  task700->add_dep(task537);
  deciq->add_task(task700);

  auto I951 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task701 = make_shared<Task701>(I845, t2, I951);
  task694->add_dep(task701);
  task701->add_dep(task537);
  deciq->add_task(task701);

  auto task702 = make_shared<Task702>(I951, f1_, t2);
  task701->add_dep(task702);
  task702->add_dep(task537);
  deciq->add_task(task702);

  auto I955 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task703 = make_shared<Task703>(I845, t2, I955);
  task694->add_dep(task703);
  task703->add_dep(task537);
  deciq->add_task(task703);

  auto task704 = make_shared<Task704>(I955, f1_, t2);
  task703->add_dep(task704);
  task704->add_dep(task537);
  deciq->add_task(task704);

  auto I959 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task705 = make_shared<Task705>(I845, t2, I959);
  task694->add_dep(task705);
  task705->add_dep(task537);
  deciq->add_task(task705);

  auto task706 = make_shared<Task706>(I959, f1_, t2);
  task705->add_dep(task706);
  task706->add_dep(task537);
  deciq->add_task(task706);

  auto I963 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task707 = make_shared<Task707>(I845, t2, I963);
  task694->add_dep(task707);
  task707->add_dep(task537);
  deciq->add_task(task707);

  auto task708 = make_shared<Task708>(I963, f1_, t2);
  task707->add_dep(task708);
  task708->add_dep(task537);
  deciq->add_task(task708);

  auto I973 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task709 = make_shared<Task709>(I845, f1_, I973);
  task694->add_dep(task709);
  task709->add_dep(task537);
  deciq->add_task(task709);

  auto task710 = make_shared<Task710>(I973, t2);
  task709->add_dep(task710);
  task710->add_dep(task537);
  deciq->add_task(task710);

  auto task711 = make_shared<Task711>(I973, t2);
  task709->add_dep(task711);
  task711->add_dep(task537);
  deciq->add_task(task711);

  auto I1005 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task712 = make_shared<Task712>(I845, f1_, I1005);
  task694->add_dep(task712);
  task712->add_dep(task537);
  deciq->add_task(task712);

  auto task713 = make_shared<Task713>(I1005, t2);
  task712->add_dep(task713);
  task713->add_dep(task537);
  deciq->add_task(task713);

  auto task714 = make_shared<Task714>(I1005, t2);
  task712->add_dep(task714);
  task714->add_dep(task537);
  deciq->add_task(task714);

  auto I1019 = make_shared<TATensor<double,4>>({virt_, virt_, active_, closed_});
  auto task715 = make_shared<Task715>(I845, t2, I1019);
  task694->add_dep(task715);
  task715->add_dep(task537);
  deciq->add_task(task715);

  auto task716 = make_shared<Task716>(I1019, f1_, t2);
  task715->add_dep(task716);
  task716->add_dep(task537);
  deciq->add_task(task716);

  auto I1023 = make_shared<TATensor<double,4>>({virt_, virt_, active_, closed_});
  auto task717 = make_shared<Task717>(I845, t2, I1023);
  task694->add_dep(task717);
  task717->add_dep(task537);
  deciq->add_task(task717);

  auto task718 = make_shared<Task718>(I1023, f1_, t2);
  task717->add_dep(task718);
  task718->add_dep(task537);
  deciq->add_task(task718);

  auto I1027 = make_shared<TATensor<double,4>>({virt_, closed_, active_, virt_});
  auto task719 = make_shared<Task719>(I845, t2, I1027);
  task694->add_dep(task719);
  task719->add_dep(task537);
  deciq->add_task(task719);

  auto task720 = make_shared<Task720>(I1027, f1_, t2);
  task719->add_dep(task720);
  task720->add_dep(task537);
  deciq->add_task(task720);

  auto I1031 = make_shared<TATensor<double,4>>({virt_, closed_, active_, virt_});
  auto task721 = make_shared<Task721>(I845, t2, I1031);
  task694->add_dep(task721);
  task721->add_dep(task537);
  deciq->add_task(task721);

  auto task722 = make_shared<Task722>(I1031, f1_, t2);
  task721->add_dep(task722);
  task722->add_dep(task537);
  deciq->add_task(task722);

  auto I1035 = make_shared<TATensor<double,4>>({closed_, virt_, active_, virt_});
  auto task723 = make_shared<Task723>(I845, t2, I1035);
  task694->add_dep(task723);
  task723->add_dep(task537);
  deciq->add_task(task723);

  auto task724 = make_shared<Task724>(I1035, f1_, t2);
  task723->add_dep(task724);
  task724->add_dep(task537);
  deciq->add_task(task724);

  auto I1039 = make_shared<TATensor<double,4>>({closed_, virt_, active_, virt_});
  auto task725 = make_shared<Task725>(I845, t2, I1039);
  task694->add_dep(task725);
  task725->add_dep(task537);
  deciq->add_task(task725);

  auto task726 = make_shared<Task726>(I1039, f1_, t2);
  task725->add_dep(task726);
  task726->add_dep(task537);
  deciq->add_task(task726);

  auto task727 = make_shared<Task727>(I845, t2, this->e0_);
  task694->add_dep(task727);
  task727->add_dep(task537);
  deciq->add_task(task727);

  auto task728 = make_shared<Task728>(I845, t2, this->e0_);
  task694->add_dep(task728);
  task728->add_dep(task537);
  deciq->add_task(task728);

  auto task729 = make_shared<Task729>(I845, v2_, t2);
  task694->add_dep(task729);
  task729->add_dep(task537);
  deciq->add_task(task729);

  auto task730 = make_shared<Task730>(I845, v2_, t2);
  task694->add_dep(task730);
  task730->add_dep(task537);
  deciq->add_task(task730);

  auto task731 = make_shared<Task731>(I845, v2_, t2);
  task694->add_dep(task731);
  task731->add_dep(task537);
  deciq->add_task(task731);

  auto task732 = make_shared<Task732>(I845, v2_, t2);
  task694->add_dep(task732);
  task732->add_dep(task537);
  deciq->add_task(task732);

  auto I1209 = make_shared<TATensor<double,4>>({active_, closed_, virt_, active_});
  auto task733 = make_shared<Task733>(I845, h1_, I1209);
  task694->add_dep(task733);
  task733->add_dep(task537);
  deciq->add_task(task733);

  auto task734 = make_shared<Task734>(I1209, t2);
  task733->add_dep(task734);
  task734->add_dep(task537);
  deciq->add_task(task734);

  auto I1212 = make_shared<TATensor<double,4>>({active_, active_, virt_, closed_});
  auto task735 = make_shared<Task735>(I845, h1_, I1212);
  task694->add_dep(task735);
  task735->add_dep(task537);
  deciq->add_task(task735);

  auto task736 = make_shared<Task736>(I1212, t2);
  task735->add_dep(task736);
  task736->add_dep(task537);
  deciq->add_task(task736);

  auto I895 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task737 = make_shared<Task737>(I698, Gamma299_(), I895);
  task538->add_dep(task737);
  task737->add_dep(task537);
  deciq->add_task(task737);

  auto I896 = make_shared<TATensor<double,4>>({active_, active_, virt_, active_});
  auto task738 = make_shared<Task738>(I895, t2, I896);
  task737->add_dep(task738);
  task738->add_dep(task537);
  deciq->add_task(task738);

  auto task739 = make_shared<Task739>(I896, f1_, t2);
  task738->add_dep(task739);
  task739->add_dep(task537);
  deciq->add_task(task739);

  auto task740 = make_shared<Task740>(I895, v2_, t2);
  task737->add_dep(task740);
  task740->add_dep(task537);
  deciq->add_task(task740);

  auto I915 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task741 = make_shared<Task741>(I698, Gamma304_(), I915);
  task538->add_dep(task741);
  task741->add_dep(task537);
  deciq->add_task(task741);

  auto I916 = make_shared<TATensor<double,4>>({active_, active_, virt_, active_});
  auto task742 = make_shared<Task742>(I915, t2, I916);
  task741->add_dep(task742);
  task742->add_dep(task537);
  deciq->add_task(task742);

  auto task743 = make_shared<Task743>(I916, t2, f1_);
  task742->add_dep(task743);
  task743->add_dep(task537);
  deciq->add_task(task743);

  auto I919 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task744 = make_shared<Task744>(I698, Gamma305_(), I919);
  task538->add_dep(task744);
  task744->add_dep(task537);
  deciq->add_task(task744);

  auto I920 = make_shared<TATensor<double,4>>({active_, virt_, active_, active_});
  auto task745 = make_shared<Task745>(I919, t2, I920);
  task744->add_dep(task745);
  task745->add_dep(task537);
  deciq->add_task(task745);

  auto task746 = make_shared<Task746>(I920, t2, f1_);
  task745->add_dep(task746);
  task746->add_dep(task537);
  deciq->add_task(task746);

  auto task747 = make_shared<Task747>(I919, v2_, t2);
  task744->add_dep(task747);
  task747->add_dep(task537);
  deciq->add_task(task747);

  auto I923 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task748 = make_shared<Task748>(I698, Gamma306_(), I923);
  task538->add_dep(task748);
  task748->add_dep(task537);
  deciq->add_task(task748);

  auto task749 = make_shared<Task749>(I923, t2);
  task748->add_dep(task749);
  task749->add_dep(task537);
  deciq->add_task(task749);

  auto I926 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task750 = make_shared<Task750>(I698, Gamma307_(), I926);
  task538->add_dep(task750);
  task750->add_dep(task537);
  deciq->add_task(task750);

  auto I927 = make_shared<TATensor<double,4>>({active_, active_, active_, virt_});
  auto task751 = make_shared<Task751>(I926, t2, I927);
  task750->add_dep(task751);
  task751->add_dep(task537);
  deciq->add_task(task751);

  auto task752 = make_shared<Task752>(I927, f1_, t2);
  task751->add_dep(task752);
  task752->add_dep(task537);
  deciq->add_task(task752);

  auto I939 = make_shared<TATensor<double,4>>({active_, active_, virt_, active_});
  auto task753 = make_shared<Task753>(I926, t2, I939);
  task750->add_dep(task753);
  task753->add_dep(task537);
  deciq->add_task(task753);

  auto task754 = make_shared<Task754>(I939, t2, f1_);
  task753->add_dep(task754);
  task754->add_dep(task537);
  deciq->add_task(task754);

  auto I1047 = make_shared<TATensor<double,4>>({active_, virt_, active_, active_});
  auto task755 = make_shared<Task755>(I926, t2, I1047);
  task750->add_dep(task755);
  task755->add_dep(task537);
  deciq->add_task(task755);

  auto task756 = make_shared<Task756>(I1047, t2, this->e0_);
  task755->add_dep(task756);
  task756->add_dep(task537);
  deciq->add_task(task756);

  auto task757 = make_shared<Task757>(I1047, f1_, t2);
  task755->add_dep(task757);
  task757->add_dep(task537);
  deciq->add_task(task757);

  auto task758 = make_shared<Task758>(I926, v2_, t2);
  task750->add_dep(task758);
  task758->add_dep(task537);
  deciq->add_task(task758);

  auto task759 = make_shared<Task759>(I926, v2_, t2);
  task750->add_dep(task759);
  task759->add_dep(task537);
  deciq->add_task(task759);

  auto I930 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task760 = make_shared<Task760>(I698, Gamma308_(), I930);
  task538->add_dep(task760);
  task760->add_dep(task537);
  deciq->add_task(task760);

  auto I931 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task761 = make_shared<Task761>(I930, t2, I931);
  task760->add_dep(task761);
  task761->add_dep(task537);
  deciq->add_task(task761);

  auto task762 = make_shared<Task762>(I931, t2, f1_);
  task761->add_dep(task762);
  task762->add_dep(task537);
  deciq->add_task(task762);

  auto task763 = make_shared<Task763>(I931, t2, f1_);
  task761->add_dep(task763);
  task763->add_dep(task537);
  deciq->add_task(task763);

  auto I997 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task764 = make_shared<Task764>(I930, t2, I997);
  task760->add_dep(task764);
  task764->add_dep(task537);
  deciq->add_task(task764);

  auto task765 = make_shared<Task765>(I997, f1_, t2);
  task764->add_dep(task765);
  task765->add_dep(task537);
  deciq->add_task(task765);

  auto I1001 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task766 = make_shared<Task766>(I930, t2, I1001);
  task760->add_dep(task766);
  task766->add_dep(task537);
  deciq->add_task(task766);

  auto task767 = make_shared<Task767>(I1001, f1_, t2);
  task766->add_dep(task767);
  task767->add_dep(task537);
  deciq->add_task(task767);

  auto I1043 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task768 = make_shared<Task768>(I930, t2, I1043);
  task760->add_dep(task768);
  task768->add_dep(task537);
  deciq->add_task(task768);

  auto task769 = make_shared<Task769>(I1043, f1_, t2);
  task768->add_dep(task769);
  task769->add_dep(task537);
  deciq->add_task(task769);

  auto task770 = make_shared<Task770>(I1043, f1_, t2);
  task768->add_dep(task770);
  task770->add_dep(task537);
  deciq->add_task(task770);

  auto I1051 = make_shared<TATensor<double,4>>({active_, active_, virt_, virt_});
  auto task771 = make_shared<Task771>(I930, t2, I1051);
  task760->add_dep(task771);
  task771->add_dep(task537);
  deciq->add_task(task771);

  auto task772 = make_shared<Task772>(I1051, t2, f1_);
  task771->add_dep(task772);
  task772->add_dep(task537);
  deciq->add_task(task772);

  auto task773 = make_shared<Task773>(I930, t2, this->e0_);
  task760->add_dep(task773);
  task773->add_dep(task537);
  deciq->add_task(task773);

  auto task774 = make_shared<Task774>(I930, v2_, t2);
  task760->add_dep(task774);
  task774->add_dep(task537);
  deciq->add_task(task774);

  auto task775 = make_shared<Task775>(I930, v2_, t2);
  task760->add_dep(task775);
  task775->add_dep(task537);
  deciq->add_task(task775);

  auto task776 = make_shared<Task776>(I930, h1_, t2);
  task760->add_dep(task776);
  task776->add_dep(task537);
  deciq->add_task(task776);

  auto task777 = make_shared<Task777>(I930, h1_, t2);
  task760->add_dep(task777);
  task777->add_dep(task537);
  deciq->add_task(task777);

  auto I966 = make_shared<TATensor<double,0>>({});
  auto task778 = make_shared<Task778>(I698, Gamma317_(), I966);
  task538->add_dep(task778);
  task778->add_dep(task537);
  deciq->add_task(task778);

  auto task779 = make_shared<Task779>(I966, t2);
  task778->add_dep(task779);
  task779->add_dep(task537);
  deciq->add_task(task779);

  auto task780 = make_shared<Task780>(I966, t2);
  task778->add_dep(task780);
  task780->add_dep(task537);
  deciq->add_task(task780);

  auto I1012 = make_shared<TATensor<double,2>>({active_, active_});
  auto task781 = make_shared<Task781>(I698, Gamma329_(), I1012);
  task538->add_dep(task781);
  task781->add_dep(task537);
  deciq->add_task(task781);

  auto task782 = make_shared<Task782>(I1012, t2);
  task781->add_dep(task782);
  task782->add_dep(task537);
  deciq->add_task(task782);

  auto task783 = make_shared<Task783>(I1012, t2);
  task781->add_dep(task783);
  task783->add_dep(task537);
  deciq->add_task(task783);

  auto I1054 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task784 = make_shared<Task784>(I698, Gamma340_(), I1054);
  task538->add_dep(task784);
  task784->add_dep(task537);
  deciq->add_task(task784);

  auto task785 = make_shared<Task785>(I1054, t2);
  task784->add_dep(task785);
  task785->add_dep(task537);
  deciq->add_task(task785);

  auto I1100 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task786 = make_shared<Task786>(I698, Gamma355_(), I1100);
  task538->add_dep(task786);
  task786->add_dep(task537);
  deciq->add_task(task786);

  auto task787 = make_shared<Task787>(I1100, v2_, t2);
  task786->add_dep(task787);
  task787->add_dep(task537);
  deciq->add_task(task787);

  auto I1154 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task788 = make_shared<Task788>(I698, Gamma373_(), I1154);
  task538->add_dep(task788);
  task788->add_dep(task537);
  deciq->add_task(task788);

  auto task789 = make_shared<Task789>(I1154, v2_, t2);
  task788->add_dep(task789);
  task789->add_dep(task537);
  deciq->add_task(task789);

  return deciq;
}


#endif
