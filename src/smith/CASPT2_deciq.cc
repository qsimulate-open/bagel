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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_deciq(const bool reset, const bool diagonal) {

  auto deciq = make_shared<Queue>();
  auto task553 = make_shared<Task553>(deci, reset);
  deciq->add_task(task553);

  auto I768 = make_shared<TATensor<double,1>>(std::vector<IndexRange>{ci_});
  auto task554 = make_shared<Task554>(deci, I768);
  task554->add_dep(task553);
  deciq->add_task(task554);

  auto I769 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task555 = make_shared<Task555>(I768, Gamma270_(), I769);
  task554->add_dep(task555);
  task555->add_dep(task553);
  deciq->add_task(task555);

  auto task556 = make_shared<Task556>(I769, t2);
  task555->add_dep(task556);
  task556->add_dep(task553);
  deciq->add_task(task556);

  auto I772 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task557 = make_shared<Task557>(I768, Gamma271_(), I772);
  task554->add_dep(task557);
  task557->add_dep(task553);
  deciq->add_task(task557);

  auto I773 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, closed_});
  auto task558 = make_shared<Task558>(I772, t2, I773);
  task557->add_dep(task558);
  task558->add_dep(task553);
  deciq->add_task(task558);

  auto task559 = make_shared<Task559>(I773, f1_, t2);
  task558->add_dep(task559);
  task559->add_dep(task553);
  deciq->add_task(task559);

  auto task560 = make_shared<Task560>(I772, t2, this->e0_);
  task557->add_dep(task560);
  task560->add_dep(task553);
  deciq->add_task(task560);

  auto task561 = make_shared<Task561>(I772, v2_, t2);
  task557->add_dep(task561);
  task561->add_dep(task553);
  deciq->add_task(task561);

  auto task562 = make_shared<Task562>(I772, v2_, t2);
  task557->add_dep(task562);
  task562->add_dep(task553);
  deciq->add_task(task562);

  auto I776 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task563 = make_shared<Task563>(I768, Gamma272_(), I776);
  task554->add_dep(task563);
  task563->add_dep(task553);
  deciq->add_task(task563);

  auto I777 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, active_});
  auto task564 = make_shared<Task564>(I776, t2, I777);
  task563->add_dep(task564);
  task564->add_dep(task553);
  deciq->add_task(task564);

  auto task565 = make_shared<Task565>(I777, f1_, t2);
  task564->add_dep(task565);
  task565->add_dep(task553);
  deciq->add_task(task565);

  auto I780 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task566 = make_shared<Task566>(I768, Gamma273_(), I780);
  task554->add_dep(task566);
  task566->add_dep(task553);
  deciq->add_task(task566);

  auto I781 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, closed_, active_});
  auto task567 = make_shared<Task567>(I780, t2, I781);
  task566->add_dep(task567);
  task567->add_dep(task553);
  deciq->add_task(task567);

  auto task568 = make_shared<Task568>(I781, t2, f1_);
  task567->add_dep(task568);
  task568->add_dep(task553);
  deciq->add_task(task568);

  auto I812 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, closed_, active_});
  auto task569 = make_shared<Task569>(I780, t2, I812);
  task566->add_dep(task569);
  task569->add_dep(task553);
  deciq->add_task(task569);

  auto task570 = make_shared<Task570>(I812, f1_, t2);
  task569->add_dep(task570);
  task570->add_dep(task553);
  deciq->add_task(task570);

  auto I784 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task571 = make_shared<Task571>(I768, Gamma274_(), I784);
  task554->add_dep(task571);
  task571->add_dep(task553);
  deciq->add_task(task571);

  auto I785 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, active_, active_});
  auto task572 = make_shared<Task572>(I784, t2, I785);
  task571->add_dep(task572);
  task572->add_dep(task553);
  deciq->add_task(task572);

  auto task573 = make_shared<Task573>(I785, t2, f1_);
  task572->add_dep(task573);
  task573->add_dep(task553);
  deciq->add_task(task573);

  auto I788 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task574 = make_shared<Task574>(I768, Gamma275_(), I788);
  task554->add_dep(task574);
  task574->add_dep(task553);
  deciq->add_task(task574);

  auto task575 = make_shared<Task575>(I788, t2);
  task574->add_dep(task575);
  task575->add_dep(task553);
  deciq->add_task(task575);

  auto I791 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task576 = make_shared<Task576>(I768, Gamma276_(), I791);
  task554->add_dep(task576);
  task576->add_dep(task553);
  deciq->add_task(task576);

  auto I792 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, closed_});
  auto task577 = make_shared<Task577>(I791, t2, I792);
  task576->add_dep(task577);
  task577->add_dep(task553);
  deciq->add_task(task577);

  auto task578 = make_shared<Task578>(I792, f1_, t2);
  task577->add_dep(task578);
  task578->add_dep(task553);
  deciq->add_task(task578);

  auto I808 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, active_, active_});
  auto task579 = make_shared<Task579>(I791, t2, I808);
  task576->add_dep(task579);
  task579->add_dep(task553);
  deciq->add_task(task579);

  auto task580 = make_shared<Task580>(I808, t2, f1_);
  task579->add_dep(task580);
  task580->add_dep(task553);
  deciq->add_task(task580);

  auto I932 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, active_});
  auto task581 = make_shared<Task581>(I791, t2, I932);
  task576->add_dep(task581);
  task581->add_dep(task553);
  deciq->add_task(task581);

  auto task582 = make_shared<Task582>(I932, t2, this->e0_);
  task581->add_dep(task582);
  task582->add_dep(task553);
  deciq->add_task(task582);

  auto task583 = make_shared<Task583>(I932, f1_, t2);
  task581->add_dep(task583);
  task583->add_dep(task553);
  deciq->add_task(task583);

  auto task584 = make_shared<Task584>(I791, v2_, t2);
  task576->add_dep(task584);
  task584->add_dep(task553);
  deciq->add_task(task584);

  auto task585 = make_shared<Task585>(I791, v2_, t2);
  task576->add_dep(task585);
  task585->add_dep(task553);
  deciq->add_task(task585);

  auto I795 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task586 = make_shared<Task586>(I768, Gamma277_(), I795);
  task554->add_dep(task586);
  task586->add_dep(task553);
  deciq->add_task(task586);

  auto I796 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task587 = make_shared<Task587>(I795, t2, I796);
  task586->add_dep(task587);
  task587->add_dep(task553);
  deciq->add_task(task587);

  auto task588 = make_shared<Task588>(I796, t2, f1_);
  task587->add_dep(task588);
  task588->add_dep(task553);
  deciq->add_task(task588);

  auto task589 = make_shared<Task589>(I796, t2, f1_);
  task587->add_dep(task589);
  task589->add_dep(task553);
  deciq->add_task(task589);

  auto I886 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, active_});
  auto task590 = make_shared<Task590>(I795, t2, I886);
  task586->add_dep(task590);
  task590->add_dep(task553);
  deciq->add_task(task590);

  auto task591 = make_shared<Task591>(I886, t2, f1_);
  task590->add_dep(task591);
  task591->add_dep(task553);
  deciq->add_task(task591);

  auto I936 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task592 = make_shared<Task592>(I795, t2, I936);
  task586->add_dep(task592);
  task592->add_dep(task553);
  deciq->add_task(task592);

  auto task593 = make_shared<Task593>(I936, t2, f1_);
  task592->add_dep(task593);
  task593->add_dep(task553);
  deciq->add_task(task593);

  auto task594 = make_shared<Task594>(I936, t2, f1_);
  task592->add_dep(task594);
  task594->add_dep(task553);
  deciq->add_task(task594);

  auto task595 = make_shared<Task595>(I795, v2_, t2);
  task586->add_dep(task595);
  task595->add_dep(task553);
  deciq->add_task(task595);

  auto task596 = make_shared<Task596>(I795, h1_, t2);
  task586->add_dep(task596);
  task596->add_dep(task553);
  deciq->add_task(task596);

  auto I803 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task597 = make_shared<Task597>(I768, Gamma279_(), I803);
  task554->add_dep(task597);
  task597->add_dep(task553);
  deciq->add_task(task597);

  auto I804 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, active_});
  auto task598 = make_shared<Task598>(I803, t2, I804);
  task597->add_dep(task598);
  task598->add_dep(task553);
  deciq->add_task(task598);

  auto task599 = make_shared<Task599>(I804, t2, f1_);
  task598->add_dep(task599);
  task599->add_dep(task553);
  deciq->add_task(task599);

  auto I815 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task600 = make_shared<Task600>(I768, Gamma282_(), I815);
  task554->add_dep(task600);
  task600->add_dep(task553);
  deciq->add_task(task600);

  auto I816 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task601 = make_shared<Task601>(I815, t2, I816);
  task600->add_dep(task601);
  task601->add_dep(task553);
  deciq->add_task(task601);

  auto task602 = make_shared<Task602>(I816, f1_, t2);
  task601->add_dep(task602);
  task602->add_dep(task553);
  deciq->add_task(task602);

  auto I820 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task603 = make_shared<Task603>(I815, t2, I820);
  task600->add_dep(task603);
  task603->add_dep(task553);
  deciq->add_task(task603);

  auto task604 = make_shared<Task604>(I820, f1_, t2);
  task603->add_dep(task604);
  task604->add_dep(task553);
  deciq->add_task(task604);

  auto I858 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task605 = make_shared<Task605>(I815, t2, I858);
  task600->add_dep(task605);
  task605->add_dep(task553);
  deciq->add_task(task605);

  auto task606 = make_shared<Task606>(I858, f1_, t2);
  task605->add_dep(task606);
  task606->add_dep(task553);
  deciq->add_task(task606);

  auto I862 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, active_});
  auto task607 = make_shared<Task607>(I815, t2, I862);
  task600->add_dep(task607);
  task607->add_dep(task553);
  deciq->add_task(task607);

  auto task608 = make_shared<Task608>(I862, f1_, t2);
  task607->add_dep(task608);
  task608->add_dep(task553);
  deciq->add_task(task608);

  auto I866 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task609 = make_shared<Task609>(I815, t2, I866);
  task600->add_dep(task609);
  task609->add_dep(task553);
  deciq->add_task(task609);

  auto task610 = make_shared<Task610>(I866, f1_, t2);
  task609->add_dep(task610);
  task610->add_dep(task553);
  deciq->add_task(task610);

  auto task611 = make_shared<Task611>(I815, v2_, t2);
  task600->add_dep(task611);
  task611->add_dep(task553);
  deciq->add_task(task611);

  auto task612 = make_shared<Task612>(I815, h1_, t2);
  task600->add_dep(task612);
  task612->add_dep(task553);
  deciq->add_task(task612);

  auto I823 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task613 = make_shared<Task613>(I768, Gamma284_(), I823);
  task554->add_dep(task613);
  task613->add_dep(task553);
  deciq->add_task(task613);

  auto task614 = make_shared<Task614>(I823, t2);
  task613->add_dep(task614);
  task614->add_dep(task553);
  deciq->add_task(task614);

  auto task615 = make_shared<Task615>(I823, t2);
  task613->add_dep(task615);
  task615->add_dep(task553);
  deciq->add_task(task615);

  auto I829 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task616 = make_shared<Task616>(I768, Gamma286_(), I829);
  task554->add_dep(task616);
  task616->add_dep(task553);
  deciq->add_task(task616);

  auto I830 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, closed_});
  auto task617 = make_shared<Task617>(I829, t2, I830);
  task616->add_dep(task617);
  task617->add_dep(task553);
  deciq->add_task(task617);

  auto task618 = make_shared<Task618>(I830, f1_, t2);
  task617->add_dep(task618);
  task618->add_dep(task553);
  deciq->add_task(task618);

  auto I834 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, closed_});
  auto task619 = make_shared<Task619>(I829, t2, I834);
  task616->add_dep(task619);
  task619->add_dep(task553);
  deciq->add_task(task619);

  auto task620 = make_shared<Task620>(I834, f1_, t2);
  task619->add_dep(task620);
  task620->add_dep(task553);
  deciq->add_task(task620);

  auto I838 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, closed_});
  auto task621 = make_shared<Task621>(I829, t2, I838);
  task616->add_dep(task621);
  task621->add_dep(task553);
  deciq->add_task(task621);

  auto task622 = make_shared<Task622>(I838, f1_, t2);
  task621->add_dep(task622);
  task622->add_dep(task553);
  deciq->add_task(task622);

  auto I842 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, closed_, virt_});
  auto task623 = make_shared<Task623>(I829, t2, I842);
  task616->add_dep(task623);
  task623->add_dep(task553);
  deciq->add_task(task623);

  auto task624 = make_shared<Task624>(I842, f1_, t2);
  task623->add_dep(task624);
  task624->add_dep(task553);
  deciq->add_task(task624);

  auto I846 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, closed_});
  auto task625 = make_shared<Task625>(I829, t2, I846);
  task616->add_dep(task625);
  task625->add_dep(task553);
  deciq->add_task(task625);

  auto task626 = make_shared<Task626>(I846, f1_, t2);
  task625->add_dep(task626);
  task626->add_dep(task553);
  deciq->add_task(task626);

  auto I850 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, closed_, virt_});
  auto task627 = make_shared<Task627>(I829, t2, I850);
  task616->add_dep(task627);
  task627->add_dep(task553);
  deciq->add_task(task627);

  auto task628 = make_shared<Task628>(I850, f1_, t2);
  task627->add_dep(task628);
  task628->add_dep(task553);
  deciq->add_task(task628);

  auto I870 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task629 = make_shared<Task629>(I829, f1_, I870);
  task616->add_dep(task629);
  task629->add_dep(task553);
  deciq->add_task(task629);

  auto task630 = make_shared<Task630>(I870, t2);
  task629->add_dep(task630);
  task630->add_dep(task553);
  deciq->add_task(task630);

  auto task631 = make_shared<Task631>(I870, t2);
  task629->add_dep(task631);
  task631->add_dep(task553);
  deciq->add_task(task631);

  auto I1013 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task632 = make_shared<Task632>(I829, f1_, I1013);
  task616->add_dep(task632);
  task632->add_dep(task553);
  deciq->add_task(task632);

  auto task633 = make_shared<Task633>(I1013, t2);
  task632->add_dep(task633);
  task633->add_dep(task553);
  deciq->add_task(task633);

  auto I1017 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task634 = make_shared<Task634>(I829, f1_, I1017);
  task616->add_dep(task634);
  task634->add_dep(task553);
  deciq->add_task(task634);

  auto task635 = make_shared<Task635>(I1017, t2);
  task634->add_dep(task635);
  task635->add_dep(task553);
  deciq->add_task(task635);

  auto task636 = make_shared<Task636>(I829, t2, this->e0_);
  task616->add_dep(task636);
  task636->add_dep(task553);
  deciq->add_task(task636);

  auto task637 = make_shared<Task637>(I829, t2, this->e0_);
  task616->add_dep(task637);
  task637->add_dep(task553);
  deciq->add_task(task637);

  auto task638 = make_shared<Task638>(I829, v2_, t2);
  task616->add_dep(task638);
  task638->add_dep(task553);
  deciq->add_task(task638);

  auto task639 = make_shared<Task639>(I829, v2_, t2);
  task616->add_dep(task639);
  task639->add_dep(task553);
  deciq->add_task(task639);

  auto task640 = make_shared<Task640>(I829, v2_, t2);
  task616->add_dep(task640);
  task640->add_dep(task553);
  deciq->add_task(task640);

  auto task641 = make_shared<Task641>(I829, v2_, t2);
  task616->add_dep(task641);
  task641->add_dep(task553);
  deciq->add_task(task641);

  auto I853 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task642 = make_shared<Task642>(I768, Gamma292_(), I853);
  task554->add_dep(task642);
  task642->add_dep(task553);
  deciq->add_task(task642);

  auto I854 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, active_});
  auto task643 = make_shared<Task643>(I853, t2, I854);
  task642->add_dep(task643);
  task643->add_dep(task553);
  deciq->add_task(task643);

  auto task644 = make_shared<Task644>(I854, f1_, t2);
  task643->add_dep(task644);
  task644->add_dep(task553);
  deciq->add_task(task644);

  auto task645 = make_shared<Task645>(I853, v2_, t2);
  task642->add_dep(task645);
  task645->add_dep(task553);
  deciq->add_task(task645);

  auto I877 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task646 = make_shared<Task646>(I768, Gamma298_(), I877);
  task554->add_dep(task646);
  task646->add_dep(task553);
  deciq->add_task(task646);

  auto I878 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, active_, active_});
  auto task647 = make_shared<Task647>(I877, t2, I878);
  task646->add_dep(task647);
  task647->add_dep(task553);
  deciq->add_task(task647);

  auto task648 = make_shared<Task648>(I878, f1_, t2);
  task647->add_dep(task648);
  task648->add_dep(task553);
  deciq->add_task(task648);

  auto I881 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task649 = make_shared<Task649>(I768, Gamma299_(), I881);
  task554->add_dep(task649);
  task649->add_dep(task553);
  deciq->add_task(task649);

  auto I882 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task650 = make_shared<Task650>(I881, t2, I882);
  task649->add_dep(task650);
  task650->add_dep(task553);
  deciq->add_task(task650);

  auto task651 = make_shared<Task651>(I882, t2, f1_);
  task650->add_dep(task651);
  task651->add_dep(task553);
  deciq->add_task(task651);

  auto task652 = make_shared<Task652>(I881, v2_, t2);
  task649->add_dep(task652);
  task652->add_dep(task553);
  deciq->add_task(task652);

  auto I889 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task653 = make_shared<Task653>(I768, Gamma301_(), I889);
  task554->add_dep(task653);
  task653->add_dep(task553);
  deciq->add_task(task653);

  auto task654 = make_shared<Task654>(I889, t2);
  task653->add_dep(task654);
  task654->add_dep(task553);
  deciq->add_task(task654);

  auto I892 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task655 = make_shared<Task655>(I768, Gamma302_(), I892);
  task554->add_dep(task655);
  task655->add_dep(task553);
  deciq->add_task(task655);

  auto I893 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, active_, closed_});
  auto task656 = make_shared<Task656>(I892, t2, I893);
  task655->add_dep(task656);
  task656->add_dep(task553);
  deciq->add_task(task656);

  auto task657 = make_shared<Task657>(I893, f1_, t2);
  task656->add_dep(task657);
  task657->add_dep(task553);
  deciq->add_task(task657);

  auto I897 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, active_, virt_});
  auto task658 = make_shared<Task658>(I892, t2, I897);
  task655->add_dep(task658);
  task658->add_dep(task553);
  deciq->add_task(task658);

  auto task659 = make_shared<Task659>(I897, f1_, t2);
  task658->add_dep(task659);
  task659->add_dep(task553);
  deciq->add_task(task659);

  auto I928 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, closed_});
  auto task660 = make_shared<Task660>(I892, t2, I928);
  task655->add_dep(task660);
  task660->add_dep(task553);
  deciq->add_task(task660);

  auto task661 = make_shared<Task661>(I928, t2, f1_);
  task660->add_dep(task661);
  task661->add_dep(task553);
  deciq->add_task(task661);

  auto I1055 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task662 = make_shared<Task662>(I892, t2, I1055);
  task655->add_dep(task662);
  task662->add_dep(task553);
  deciq->add_task(task662);

  auto task663 = make_shared<Task663>(I1055, t2, this->e0_);
  task662->add_dep(task663);
  task663->add_dep(task553);
  deciq->add_task(task663);

  auto task664 = make_shared<Task664>(I1055, f1_, t2);
  task662->add_dep(task664);
  task664->add_dep(task553);
  deciq->add_task(task664);

  auto task665 = make_shared<Task665>(I892, v2_, t2);
  task655->add_dep(task665);
  task665->add_dep(task553);
  deciq->add_task(task665);

  auto task666 = make_shared<Task666>(I892, v2_, t2);
  task655->add_dep(task666);
  task666->add_dep(task553);
  deciq->add_task(task666);

  auto I900 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task667 = make_shared<Task667>(I768, Gamma304_(), I900);
  task554->add_dep(task667);
  task667->add_dep(task553);
  deciq->add_task(task667);

  auto task668 = make_shared<Task668>(I900, t2);
  task667->add_dep(task668);
  task668->add_dep(task553);
  deciq->add_task(task668);

  auto task669 = make_shared<Task669>(I900, t2);
  task667->add_dep(task669);
  task669->add_dep(task553);
  deciq->add_task(task669);

  auto task670 = make_shared<Task670>(I900, t2);
  task667->add_dep(task670);
  task670->add_dep(task553);
  deciq->add_task(task670);

  auto I903 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task671 = make_shared<Task671>(I768, Gamma305_(), I903);
  task554->add_dep(task671);
  task671->add_dep(task553);
  deciq->add_task(task671);

  auto I904 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, active_, closed_});
  auto task672 = make_shared<Task672>(I903, t2, I904);
  task671->add_dep(task672);
  task672->add_dep(task553);
  deciq->add_task(task672);

  auto task673 = make_shared<Task673>(I904, f1_, t2);
  task672->add_dep(task673);
  task673->add_dep(task553);
  deciq->add_task(task673);

  auto I908 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, active_, virt_});
  auto task674 = make_shared<Task674>(I903, t2, I908);
  task671->add_dep(task674);
  task674->add_dep(task553);
  deciq->add_task(task674);

  auto task675 = make_shared<Task675>(I908, f1_, t2);
  task674->add_dep(task675);
  task675->add_dep(task553);
  deciq->add_task(task675);

  auto task676 = make_shared<Task676>(I908, f1_, t2);
  task674->add_dep(task676);
  task676->add_dep(task553);
  deciq->add_task(task676);

  auto I924 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, virt_});
  auto task677 = make_shared<Task677>(I903, t2, I924);
  task671->add_dep(task677);
  task677->add_dep(task553);
  deciq->add_task(task677);

  auto task678 = make_shared<Task678>(I924, t2, f1_);
  task677->add_dep(task678);
  task678->add_dep(task553);
  deciq->add_task(task678);

  auto I947 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, closed_});
  auto task679 = make_shared<Task679>(I903, t2, I947);
  task671->add_dep(task679);
  task679->add_dep(task553);
  deciq->add_task(task679);

  auto task680 = make_shared<Task680>(I947, f1_, t2);
  task679->add_dep(task680);
  task680->add_dep(task553);
  deciq->add_task(task680);

  auto I951 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, virt_});
  auto task681 = make_shared<Task681>(I903, t2, I951);
  task671->add_dep(task681);
  task681->add_dep(task553);
  deciq->add_task(task681);

  auto task682 = make_shared<Task682>(I951, f1_, t2);
  task681->add_dep(task682);
  task682->add_dep(task553);
  deciq->add_task(task682);

  auto I958 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, closed_});
  auto task683 = make_shared<Task683>(I903, t2, I958);
  task671->add_dep(task683);
  task683->add_dep(task553);
  deciq->add_task(task683);

  auto task684 = make_shared<Task684>(I958, f1_, t2);
  task683->add_dep(task684);
  task684->add_dep(task553);
  deciq->add_task(task684);

  auto I962 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, virt_});
  auto task685 = make_shared<Task685>(I903, t2, I962);
  task671->add_dep(task685);
  task685->add_dep(task553);
  deciq->add_task(task685);

  auto task686 = make_shared<Task686>(I962, f1_, t2);
  task685->add_dep(task686);
  task686->add_dep(task553);
  deciq->add_task(task686);

  auto I978 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, closed_, virt_});
  auto task687 = make_shared<Task687>(I903, t2, I978);
  task671->add_dep(task687);
  task687->add_dep(task553);
  deciq->add_task(task687);

  auto task688 = make_shared<Task688>(I978, t2, f1_);
  task687->add_dep(task688);
  task688->add_dep(task553);
  deciq->add_task(task688);

  auto task689 = make_shared<Task689>(I978, t2, f1_);
  task687->add_dep(task689);
  task689->add_dep(task553);
  deciq->add_task(task689);

  auto I1051 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task690 = make_shared<Task690>(I903, t2, I1051);
  task671->add_dep(task690);
  task690->add_dep(task553);
  deciq->add_task(task690);

  auto task691 = make_shared<Task691>(I1051, f1_, t2);
  task690->add_dep(task691);
  task691->add_dep(task553);
  deciq->add_task(task691);

  auto I1063 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task692 = make_shared<Task692>(I903, t2, I1063);
  task671->add_dep(task692);
  task692->add_dep(task553);
  deciq->add_task(task692);

  auto task693 = make_shared<Task693>(I1063, t2, this->e0_);
  task692->add_dep(task693);
  task693->add_dep(task553);
  deciq->add_task(task693);

  auto task694 = make_shared<Task694>(I1063, f1_, t2);
  task692->add_dep(task694);
  task694->add_dep(task553);
  deciq->add_task(task694);

  auto task695 = make_shared<Task695>(I903, t2, this->e0_);
  task671->add_dep(task695);
  task695->add_dep(task553);
  deciq->add_task(task695);

  auto task696 = make_shared<Task696>(I903, t2, this->e0_);
  task671->add_dep(task696);
  task696->add_dep(task553);
  deciq->add_task(task696);

  auto task697 = make_shared<Task697>(I903, v2_, t2);
  task671->add_dep(task697);
  task697->add_dep(task553);
  deciq->add_task(task697);

  auto task698 = make_shared<Task698>(I903, v2_, t2);
  task671->add_dep(task698);
  task698->add_dep(task553);
  deciq->add_task(task698);

  auto task699 = make_shared<Task699>(I903, v2_, t2);
  task671->add_dep(task699);
  task699->add_dep(task553);
  deciq->add_task(task699);

  auto task700 = make_shared<Task700>(I903, v2_, t2);
  task671->add_dep(task700);
  task700->add_dep(task553);
  deciq->add_task(task700);

  auto task701 = make_shared<Task701>(I903, v2_, t2);
  task671->add_dep(task701);
  task701->add_dep(task553);
  deciq->add_task(task701);

  auto task702 = make_shared<Task702>(I903, v2_, t2);
  task671->add_dep(task702);
  task702->add_dep(task553);
  deciq->add_task(task702);

  auto task703 = make_shared<Task703>(I903, v2_, t2);
  task671->add_dep(task703);
  task703->add_dep(task553);
  deciq->add_task(task703);

  auto task704 = make_shared<Task704>(I903, v2_, t2);
  task671->add_dep(task704);
  task704->add_dep(task553);
  deciq->add_task(task704);

  auto task705 = make_shared<Task705>(I903, v2_, t2);
  task671->add_dep(task705);
  task705->add_dep(task553);
  deciq->add_task(task705);

  auto task706 = make_shared<Task706>(I903, v2_, t2);
  task671->add_dep(task706);
  task706->add_dep(task553);
  deciq->add_task(task706);

  auto I911 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task707 = make_shared<Task707>(I768, Gamma307_(), I911);
  task554->add_dep(task707);
  task707->add_dep(task553);
  deciq->add_task(task707);

  auto I912 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, active_, active_});
  auto task708 = make_shared<Task708>(I911, t2, I912);
  task707->add_dep(task708);
  task708->add_dep(task553);
  deciq->add_task(task708);

  auto task709 = make_shared<Task709>(I912, f1_, t2);
  task708->add_dep(task709);
  task709->add_dep(task553);
  deciq->add_task(task709);

  auto I915 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task710 = make_shared<Task710>(I768, Gamma308_(), I915);
  task554->add_dep(task710);
  task710->add_dep(task553);
  deciq->add_task(task710);

  auto I916 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task711 = make_shared<Task711>(I915, t2, I916);
  task710->add_dep(task711);
  task711->add_dep(task553);
  deciq->add_task(task711);

  auto task712 = make_shared<Task712>(I916, t2, f1_);
  task711->add_dep(task712);
  task712->add_dep(task553);
  deciq->add_task(task712);

  auto task713 = make_shared<Task713>(I916, t2, f1_);
  task711->add_dep(task713);
  task713->add_dep(task553);
  deciq->add_task(task713);

  auto I970 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task714 = make_shared<Task714>(I915, t2, I970);
  task710->add_dep(task714);
  task714->add_dep(task553);
  deciq->add_task(task714);

  auto task715 = make_shared<Task715>(I970, t2, f1_);
  task714->add_dep(task715);
  task715->add_dep(task553);
  deciq->add_task(task715);

  auto task716 = make_shared<Task716>(I970, t2, f1_);
  task714->add_dep(task716);
  task716->add_dep(task553);
  deciq->add_task(task716);

  auto I1021 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task717 = make_shared<Task717>(I915, t2, I1021);
  task710->add_dep(task717);
  task717->add_dep(task553);
  deciq->add_task(task717);

  auto task718 = make_shared<Task718>(I1021, f1_, t2);
  task717->add_dep(task718);
  task718->add_dep(task553);
  deciq->add_task(task718);

  auto I1025 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task719 = make_shared<Task719>(I915, t2, I1025);
  task710->add_dep(task719);
  task719->add_dep(task553);
  deciq->add_task(task719);

  auto task720 = make_shared<Task720>(I1025, f1_, t2);
  task719->add_dep(task720);
  task720->add_dep(task553);
  deciq->add_task(task720);

  auto I1029 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task721 = make_shared<Task721>(I915, t2, I1029);
  task710->add_dep(task721);
  task721->add_dep(task553);
  deciq->add_task(task721);

  auto task722 = make_shared<Task722>(I1029, f1_, t2);
  task721->add_dep(task722);
  task722->add_dep(task553);
  deciq->add_task(task722);

  auto I1033 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task723 = make_shared<Task723>(I915, t2, I1033);
  task710->add_dep(task723);
  task723->add_dep(task553);
  deciq->add_task(task723);

  auto task724 = make_shared<Task724>(I1033, f1_, t2);
  task723->add_dep(task724);
  task724->add_dep(task553);
  deciq->add_task(task724);

  auto I1043 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task725 = make_shared<Task725>(I915, f1_, I1043);
  task710->add_dep(task725);
  task725->add_dep(task553);
  deciq->add_task(task725);

  auto task726 = make_shared<Task726>(I1043, t2);
  task725->add_dep(task726);
  task726->add_dep(task553);
  deciq->add_task(task726);

  auto task727 = make_shared<Task727>(I1043, t2);
  task725->add_dep(task727);
  task727->add_dep(task553);
  deciq->add_task(task727);

  auto I1075 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task728 = make_shared<Task728>(I915, f1_, I1075);
  task710->add_dep(task728);
  task728->add_dep(task553);
  deciq->add_task(task728);

  auto task729 = make_shared<Task729>(I1075, t2);
  task728->add_dep(task729);
  task729->add_dep(task553);
  deciq->add_task(task729);

  auto task730 = make_shared<Task730>(I1075, t2);
  task728->add_dep(task730);
  task730->add_dep(task553);
  deciq->add_task(task730);

  auto I1089 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, closed_});
  auto task731 = make_shared<Task731>(I915, t2, I1089);
  task710->add_dep(task731);
  task731->add_dep(task553);
  deciq->add_task(task731);

  auto task732 = make_shared<Task732>(I1089, f1_, t2);
  task731->add_dep(task732);
  task732->add_dep(task553);
  deciq->add_task(task732);

  auto I1093 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, closed_});
  auto task733 = make_shared<Task733>(I915, t2, I1093);
  task710->add_dep(task733);
  task733->add_dep(task553);
  deciq->add_task(task733);

  auto task734 = make_shared<Task734>(I1093, f1_, t2);
  task733->add_dep(task734);
  task734->add_dep(task553);
  deciq->add_task(task734);

  auto I1097 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, virt_});
  auto task735 = make_shared<Task735>(I915, t2, I1097);
  task710->add_dep(task735);
  task735->add_dep(task553);
  deciq->add_task(task735);

  auto task736 = make_shared<Task736>(I1097, f1_, t2);
  task735->add_dep(task736);
  task736->add_dep(task553);
  deciq->add_task(task736);

  auto I1101 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, virt_});
  auto task737 = make_shared<Task737>(I915, t2, I1101);
  task710->add_dep(task737);
  task737->add_dep(task553);
  deciq->add_task(task737);

  auto task738 = make_shared<Task738>(I1101, f1_, t2);
  task737->add_dep(task738);
  task738->add_dep(task553);
  deciq->add_task(task738);

  auto I1105 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, virt_});
  auto task739 = make_shared<Task739>(I915, t2, I1105);
  task710->add_dep(task739);
  task739->add_dep(task553);
  deciq->add_task(task739);

  auto task740 = make_shared<Task740>(I1105, f1_, t2);
  task739->add_dep(task740);
  task740->add_dep(task553);
  deciq->add_task(task740);

  auto I1109 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, virt_});
  auto task741 = make_shared<Task741>(I915, t2, I1109);
  task710->add_dep(task741);
  task741->add_dep(task553);
  deciq->add_task(task741);

  auto task742 = make_shared<Task742>(I1109, f1_, t2);
  task741->add_dep(task742);
  task742->add_dep(task553);
  deciq->add_task(task742);

  auto task743 = make_shared<Task743>(I915, t2, this->e0_);
  task710->add_dep(task743);
  task743->add_dep(task553);
  deciq->add_task(task743);

  auto task744 = make_shared<Task744>(I915, t2, this->e0_);
  task710->add_dep(task744);
  task744->add_dep(task553);
  deciq->add_task(task744);

  auto task745 = make_shared<Task745>(I915, v2_, t2);
  task710->add_dep(task745);
  task745->add_dep(task553);
  deciq->add_task(task745);

  auto task746 = make_shared<Task746>(I915, v2_, t2);
  task710->add_dep(task746);
  task746->add_dep(task553);
  deciq->add_task(task746);

  auto task747 = make_shared<Task747>(I915, v2_, t2);
  task710->add_dep(task747);
  task747->add_dep(task553);
  deciq->add_task(task747);

  auto task748 = make_shared<Task748>(I915, v2_, t2);
  task710->add_dep(task748);
  task748->add_dep(task553);
  deciq->add_task(task748);

  auto I1279 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, closed_, virt_, active_});
  auto task749 = make_shared<Task749>(I915, h1_, I1279);
  task710->add_dep(task749);
  task749->add_dep(task553);
  deciq->add_task(task749);

  auto task750 = make_shared<Task750>(I1279, t2);
  task749->add_dep(task750);
  task750->add_dep(task553);
  deciq->add_task(task750);

  auto I1282 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, closed_});
  auto task751 = make_shared<Task751>(I915, h1_, I1282);
  task710->add_dep(task751);
  task751->add_dep(task553);
  deciq->add_task(task751);

  auto task752 = make_shared<Task752>(I1282, t2);
  task751->add_dep(task752);
  task752->add_dep(task553);
  deciq->add_task(task752);

  auto I965 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task753 = make_shared<Task753>(I768, Gamma321_(), I965);
  task554->add_dep(task753);
  task753->add_dep(task553);
  deciq->add_task(task753);

  auto I966 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, active_});
  auto task754 = make_shared<Task754>(I965, t2, I966);
  task753->add_dep(task754);
  task754->add_dep(task553);
  deciq->add_task(task754);

  auto task755 = make_shared<Task755>(I966, f1_, t2);
  task754->add_dep(task755);
  task755->add_dep(task553);
  deciq->add_task(task755);

  auto task756 = make_shared<Task756>(I965, v2_, t2);
  task753->add_dep(task756);
  task756->add_dep(task553);
  deciq->add_task(task756);

  auto I985 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task757 = make_shared<Task757>(I768, Gamma326_(), I985);
  task554->add_dep(task757);
  task757->add_dep(task553);
  deciq->add_task(task757);

  auto I986 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, active_});
  auto task758 = make_shared<Task758>(I985, t2, I986);
  task757->add_dep(task758);
  task758->add_dep(task553);
  deciq->add_task(task758);

  auto task759 = make_shared<Task759>(I986, t2, f1_);
  task758->add_dep(task759);
  task759->add_dep(task553);
  deciq->add_task(task759);

  auto I989 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task760 = make_shared<Task760>(I768, Gamma327_(), I989);
  task554->add_dep(task760);
  task760->add_dep(task553);
  deciq->add_task(task760);

  auto I990 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, active_, active_});
  auto task761 = make_shared<Task761>(I989, t2, I990);
  task760->add_dep(task761);
  task761->add_dep(task553);
  deciq->add_task(task761);

  auto task762 = make_shared<Task762>(I990, t2, f1_);
  task761->add_dep(task762);
  task762->add_dep(task553);
  deciq->add_task(task762);

  auto task763 = make_shared<Task763>(I989, v2_, t2);
  task760->add_dep(task763);
  task763->add_dep(task553);
  deciq->add_task(task763);

  auto I993 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task764 = make_shared<Task764>(I768, Gamma328_(), I993);
  task554->add_dep(task764);
  task764->add_dep(task553);
  deciq->add_task(task764);

  auto task765 = make_shared<Task765>(I993, t2);
  task764->add_dep(task765);
  task765->add_dep(task553);
  deciq->add_task(task765);

  auto I996 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task766 = make_shared<Task766>(I768, Gamma329_(), I996);
  task554->add_dep(task766);
  task766->add_dep(task553);
  deciq->add_task(task766);

  auto I997 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, virt_});
  auto task767 = make_shared<Task767>(I996, t2, I997);
  task766->add_dep(task767);
  task767->add_dep(task553);
  deciq->add_task(task767);

  auto task768 = make_shared<Task768>(I997, f1_, t2);
  task767->add_dep(task768);
  task768->add_dep(task553);
  deciq->add_task(task768);

  auto I1009 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, active_});
  auto task769 = make_shared<Task769>(I996, t2, I1009);
  task766->add_dep(task769);
  task769->add_dep(task553);
  deciq->add_task(task769);

  auto task770 = make_shared<Task770>(I1009, t2, f1_);
  task769->add_dep(task770);
  task770->add_dep(task553);
  deciq->add_task(task770);

  auto I1117 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, active_, active_});
  auto task771 = make_shared<Task771>(I996, t2, I1117);
  task766->add_dep(task771);
  task771->add_dep(task553);
  deciq->add_task(task771);

  auto task772 = make_shared<Task772>(I1117, t2, this->e0_);
  task771->add_dep(task772);
  task772->add_dep(task553);
  deciq->add_task(task772);

  auto task773 = make_shared<Task773>(I1117, f1_, t2);
  task771->add_dep(task773);
  task773->add_dep(task553);
  deciq->add_task(task773);

  auto task774 = make_shared<Task774>(I996, v2_, t2);
  task766->add_dep(task774);
  task774->add_dep(task553);
  deciq->add_task(task774);

  auto task775 = make_shared<Task775>(I996, v2_, t2);
  task766->add_dep(task775);
  task775->add_dep(task553);
  deciq->add_task(task775);

  auto I1000 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task776 = make_shared<Task776>(I768, Gamma330_(), I1000);
  task554->add_dep(task776);
  task776->add_dep(task553);
  deciq->add_task(task776);

  auto I1001 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task777 = make_shared<Task777>(I1000, t2, I1001);
  task776->add_dep(task777);
  task777->add_dep(task553);
  deciq->add_task(task777);

  auto task778 = make_shared<Task778>(I1001, t2, f1_);
  task777->add_dep(task778);
  task778->add_dep(task553);
  deciq->add_task(task778);

  auto task779 = make_shared<Task779>(I1001, t2, f1_);
  task777->add_dep(task779);
  task779->add_dep(task553);
  deciq->add_task(task779);

  auto I1067 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task780 = make_shared<Task780>(I1000, t2, I1067);
  task776->add_dep(task780);
  task780->add_dep(task553);
  deciq->add_task(task780);

  auto task781 = make_shared<Task781>(I1067, f1_, t2);
  task780->add_dep(task781);
  task781->add_dep(task553);
  deciq->add_task(task781);

  auto I1071 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task782 = make_shared<Task782>(I1000, t2, I1071);
  task776->add_dep(task782);
  task782->add_dep(task553);
  deciq->add_task(task782);

  auto task783 = make_shared<Task783>(I1071, f1_, t2);
  task782->add_dep(task783);
  task783->add_dep(task553);
  deciq->add_task(task783);

  auto I1113 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task784 = make_shared<Task784>(I1000, t2, I1113);
  task776->add_dep(task784);
  task784->add_dep(task553);
  deciq->add_task(task784);

  auto task785 = make_shared<Task785>(I1113, f1_, t2);
  task784->add_dep(task785);
  task785->add_dep(task553);
  deciq->add_task(task785);

  auto task786 = make_shared<Task786>(I1113, f1_, t2);
  task784->add_dep(task786);
  task786->add_dep(task553);
  deciq->add_task(task786);

  auto I1121 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, virt_, virt_});
  auto task787 = make_shared<Task787>(I1000, t2, I1121);
  task776->add_dep(task787);
  task787->add_dep(task553);
  deciq->add_task(task787);

  auto task788 = make_shared<Task788>(I1121, t2, f1_);
  task787->add_dep(task788);
  task788->add_dep(task553);
  deciq->add_task(task788);

  auto task789 = make_shared<Task789>(I1000, t2, this->e0_);
  task776->add_dep(task789);
  task789->add_dep(task553);
  deciq->add_task(task789);

  auto task790 = make_shared<Task790>(I1000, v2_, t2);
  task776->add_dep(task790);
  task790->add_dep(task553);
  deciq->add_task(task790);

  auto task791 = make_shared<Task791>(I1000, v2_, t2);
  task776->add_dep(task791);
  task791->add_dep(task553);
  deciq->add_task(task791);

  auto task792 = make_shared<Task792>(I1000, h1_, t2);
  task776->add_dep(task792);
  task792->add_dep(task553);
  deciq->add_task(task792);

  auto task793 = make_shared<Task793>(I1000, h1_, t2);
  task776->add_dep(task793);
  task793->add_dep(task553);
  deciq->add_task(task793);

  auto I1036 = make_shared<TATensor<double,0>>(std::vector<IndexRange>{});
  auto task794 = make_shared<Task794>(I768, Gamma339_(), I1036);
  task554->add_dep(task794);
  task794->add_dep(task553);
  deciq->add_task(task794);

  auto task795 = make_shared<Task795>(I1036, t2);
  task794->add_dep(task795);
  task795->add_dep(task553);
  deciq->add_task(task795);

  auto task796 = make_shared<Task796>(I1036, t2);
  task794->add_dep(task796);
  task796->add_dep(task553);
  deciq->add_task(task796);

  auto I1082 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task797 = make_shared<Task797>(I768, Gamma351_(), I1082);
  task554->add_dep(task797);
  task797->add_dep(task553);
  deciq->add_task(task797);

  auto task798 = make_shared<Task798>(I1082, t2);
  task797->add_dep(task798);
  task798->add_dep(task553);
  deciq->add_task(task798);

  auto task799 = make_shared<Task799>(I1082, t2);
  task797->add_dep(task799);
  task799->add_dep(task553);
  deciq->add_task(task799);

  auto I1124 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task800 = make_shared<Task800>(I768, Gamma362_(), I1124);
  task554->add_dep(task800);
  task800->add_dep(task553);
  deciq->add_task(task800);

  auto task801 = make_shared<Task801>(I1124, t2);
  task800->add_dep(task801);
  task801->add_dep(task553);
  deciq->add_task(task801);

  auto I1170 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task802 = make_shared<Task802>(I768, Gamma377_(), I1170);
  task554->add_dep(task802);
  task802->add_dep(task553);
  deciq->add_task(task802);

  auto task803 = make_shared<Task803>(I1170, v2_, t2);
  task802->add_dep(task803);
  task803->add_dep(task553);
  deciq->add_task(task803);

  auto I1224 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task804 = make_shared<Task804>(I768, Gamma395_(), I1224);
  task554->add_dep(task804);
  task804->add_dep(task553);
  deciq->add_task(task804);

  auto task805 = make_shared<Task805>(I1224, v2_, t2);
  task804->add_dep(task805);
  task805->add_dep(task553);
  deciq->add_task(task805);

  return deciq;
}


#endif
