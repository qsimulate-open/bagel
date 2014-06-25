//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_compute_cispaces.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <tuple>

#include <src/dimer/dimer.h>

using namespace std;
using namespace bagel;

void Dimer::get_spaces(shared_ptr<const PTree> idata, vector<vector<int>>& spaces_A, vector<vector<int>>& spaces_B) {
  auto space = idata->get_child_optional("space");
  if (space) {
    for (auto& s : *space) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    spaces_B = spaces_A;
  }
  else {
    auto spacea = idata->get_child_optional("space_a");
    auto spaceb = idata->get_child_optional("space_b");
    if (!(spacea && spaceb)) {
      throw runtime_error("Must specify either space keywords or BOTH space_a and space_b");
    }
    for (auto& s : *spacea) { spaces_A.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
    for (auto& s : *spaceb) { spaces_B.push_back(vector<int>{s->get<int>("charge"), s->get<int>("spin"), s->get<int>("nstate")}); }
  }
}

shared_ptr<DimerCAS> Dimer::compute_cispace(const shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  auto d1 = make_shared<Determinants>(active_refs_.first->nact(), nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(active_refs_.second->nact(), nelea.second, neleb.second, /*compress*/false, /*mute*/true);
  auto out = make_shared<DimerCAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer castime;

  shared_ptr<const PTree> fcidata = idata->get_child_optional("fci");
  if (!fcidata) fcidata = make_shared<const PTree>();

  // Embedded CAS-CI calculations
  cout << "    Starting embedded CAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_casci<0>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded CAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Charge, spin and number of states needs to be specified for each space");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_casci<1>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerDistCAS> Dimer::compute_distcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  auto d1 = make_shared<Determinants>(active_refs_.first->nact(), nelea.first, neleb.first, /*compress*/false, /*mute*/true);
  auto d2 = make_shared<Determinants>(active_refs_.second->nact(), nelea.second, neleb.second, /*compress*/false, /*mute*/true);
  auto out = make_shared<DimerDistCAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer castime;

  shared_ptr<const PTree> fcidata = idata->get_child_optional("fci");
  if (!fcidata) fcidata = make_shared<const PTree>();

  // Embedded CAS-CI calculations
  cout << "    Starting embedded distributed CAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_distcasci<0>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded distributed CAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Charge, spin and number of states needs to be specified for each space");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_distcasci<1>(fcidata, charge, spin, nstate));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerRAS> Dimer::compute_rcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  // { {nras1, nras2, nras3}, max holes, max particles }
  pair<tuple<array<int, 3>, int, int>, tuple<array<int, 3>, int, int>> ras_desc;

  // Sample:
  // "restricted" : [ { "orbitals" : [1, 2, 3], "max_holes" : 0, "max_particles" : 2 } ],
  //  puts 1 orbital in RAS1 with no holes allowed, 2 orbital in RAS2, and 3 orbitals in RAS3 with up to 2 particles
  auto restrictions = idata->get_child("restricted");

  auto get_restricted_data = [] (shared_ptr<const PTree> i) {
    return make_tuple(i->get_array<int, 3>("orbitals"), i->get<int>("max_holes"), i->get<int>("max_particles"));
  };

  if (restrictions->size() == 1) {
    ras_desc = make_pair( get_restricted_data(*restrictions->begin()), get_restricted_data(*restrictions->begin()) );
  }
  else if (restrictions->size() == 2) {
    auto iter = restrictions->begin();
    auto tmp1 = get_restricted_data(*iter++);
    auto tmp2 = get_restricted_data(*iter);
    ras_desc = make_pair(tmp1, tmp2);
  }
  else throw logic_error("One or two sets of restrictions must be provided.");

  // This is less than ideal. It'd be better to have some sort of generator object that can be passed around.
  auto d1 = make_shared<RASDeterminants>(get<0>(ras_desc.first), nelea.first, neleb.first, get<1>(ras_desc.first), get<2>(ras_desc.first), true);
  auto d2 = make_shared<RASDeterminants>(get<0>(ras_desc.second), nelea.second, neleb.second, get<1>(ras_desc.second), get<2>(ras_desc.second), true);

  auto out = make_shared<DimerRAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer castime;

  shared_ptr<const PTree> rasdata = idata->get_child_optional("ras");
  if (!rasdata) rasdata = make_shared<const PTree>();

  // Embedded RAS-CI calculations
  cout << "    Starting embedded RAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_rasci<0>(rasdata, charge, spin, nstate, ras_desc.first));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded RAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_rasci<1>(rasdata, charge, spin, nstate, ras_desc.second));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}


shared_ptr<DimerDistRAS> Dimer::compute_distrcispace(const std::shared_ptr<const PTree> idata) {
  embed_refs();
  pair<int,int> nelea = make_pair(isolated_refs_.first->nclosed() - active_refs_.first->nclosed(),
                                  isolated_refs_.second->nclosed() - active_refs_.second->nclosed());
  pair<int,int> neleb = nelea;

  // { {nras1, nras2, nras3}, max holes, max particles }
  pair<tuple<array<int, 3>, int, int>, tuple<array<int, 3>, int, int>> ras_desc;

  // Sample:
  // "restricted" : [ { "orbitals" : [1, 2, 3], "max_holes" : 0, "max_particles" : 2 } ],
  //  puts 1 orbital in RAS1 with no holes allowed, 2 orbital in RAS2, and 3 orbitals in RAS3 with up to 2 particles
  auto restrictions = idata->get_child("restricted");

  auto get_restricted_data = [] (shared_ptr<const PTree> i) {
    return make_tuple(i->get_array<int, 3>("orbitals"), i->get<int>("max_holes"), i->get<int>("max_particles"));
  };

  if (restrictions->size() == 1) {
    ras_desc = make_pair( get_restricted_data(*restrictions->begin()), get_restricted_data(*restrictions->begin()) );
  }
  else if (restrictions->size() == 2) {
    auto iter = restrictions->begin();
    auto tmp1 = get_restricted_data(*iter++);
    auto tmp2 = get_restricted_data(*iter);
    ras_desc = make_pair(tmp1, tmp2);
  }
  else throw logic_error("One or two sets of restrictions must be provided.");

  // This is less than ideal. It'd be better to have some sort of generator object that can be passed around.
  auto d1 = make_shared<RASDeterminants>(get<0>(ras_desc.first), nelea.first, neleb.first, get<1>(ras_desc.first), get<2>(ras_desc.first), true);
  auto d2 = make_shared<RASDeterminants>(get<0>(ras_desc.second), nelea.second, neleb.second, get<1>(ras_desc.second), get<2>(ras_desc.second), true);

  auto out = make_shared<DimerDistRAS>(make_pair(d1, d2), nelea, neleb);

  vector<vector<int>> spaces_A, spaces_B;
  get_spaces(idata, spaces_A, spaces_B);

  Timer castime;

  shared_ptr<const PTree> rasdata = idata->get_child_optional("ras");
  if (!rasdata) rasdata = make_shared<const PTree>();

  // Embedded RAS-CI calculations
  cout << "    Starting embedded RAS-CI calculations on monomer A" << endl;
  for (auto& ispace : spaces_A) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<0>(embedded_distrasci<0>(rasdata, charge, spin, nstate, ras_desc.first));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }

  cout << endl << "    Starting embedded RAS-CI calculations on monomer B" << endl;
  for (auto& ispace : spaces_B) {
    if (ispace.size() != 3) throw runtime_error("Spaces should be input as \"space = charge, spin, nstates\"");
    const int charge = ispace.at(0);
    const int spin = ispace.at(1);
    const int nstate = ispace.at(2);

    out->insert<1>(embedded_distrasci<1>(rasdata, charge, spin, nstate, ras_desc.second));

    cout << "      - charge: " << charge << ", spin: " << spin << ", nstates: " << nstate
                               << fixed << setw(10) << setprecision(2) << castime.tick() << endl;
  }


  return out;
}
