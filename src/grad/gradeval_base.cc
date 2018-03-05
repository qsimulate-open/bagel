//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradeval_base.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/grad/gradeval_base.h>
#include <src/grad/dkhgrad.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/resources.h>
#include <src/util/parallel/mpi_interface.h>
#include <array>

using namespace std;
using namespace bagel;

shared_ptr<GradFile> GradEval_base::contract_gradient(const shared_ptr<const Matrix> d, const shared_ptr<const Matrix> w,
                                                      const shared_ptr<const DFDist> o, const shared_ptr<const Matrix> o2,
                                                      const shared_ptr<const Matrix> v, const bool numerical,
                                                      const shared_ptr<const Geometry> g2, const shared_ptr<const DFDist> g2o, const shared_ptr<const Matrix> g2o2) {
  grad_->zero();

  if (!numerical) {
    vector<shared_ptr<GradTask>> task  = contract_grad2e(o);

    vector<shared_ptr<GradTask>> task2;
    if (geom_->hcoreinfo()->dkh()) {
      auto dkh = make_shared<DKHgrad>(geom_);
      task2 = contract_graddkh1e(dkh->compute(d, w));
    } else {
      task2 = contract_grad1e<GradTask1>(d, w);
    }

    vector<shared_ptr<GradTask>> task3 = contract_grad2e_2index(o2);
    task.insert(task.end(), task2.begin(), task2.end());
    task.insert(task.end(), task3.begin(), task3.end());
    if (v) {
      vector<shared_ptr<GradTask>> task4 = contract_grad1e<GradTask1s>(v, v);
      task.insert(task.end(), task4.begin(), task4.end());
    }

    if (g2 && g2o) {
      vector<shared_ptr<GradTask>> task0 = contract_grad2e(g2o, g2);
      task.insert(task.end(), task0.begin(), task0.end());
    }
    if (g2 && g2o2) {
      vector<shared_ptr<GradTask>> task0 = contract_grad2e_2index(g2o2, g2);
      task.insert(task.end(), task0.begin(), task0.end());
    }

    TaskQueue<shared_ptr<GradTask>> tq(move(task));
    tq.compute();
  } else {
    vector<shared_ptr<GradTask>> task = contract_grad1e<GradTask1s>(v, v);

    TaskQueue<shared_ptr<GradTask>> tq(move(task));
    tq.compute();
  }

  if (!v)
    *grad_ += *geom_->compute_grad_vnuc();

  grad_->allreduce();

  return grad_;
}


template<typename TaskType>
vector<shared_ptr<GradTask>> GradEval_base::contract_grad1e(const shared_ptr<const Matrix> d, const shared_ptr<const Matrix> w) {
  return contract_grad1e<TaskType>(d, d, w);
}


template<typename TaskType>
vector<shared_ptr<GradTask>> GradEval_base::contract_grad1e(const shared_ptr<const Matrix> nmat, const shared_ptr<const Matrix> kmat, const shared_ptr<const Matrix> omat) {
  vector<shared_ptr<GradTask>> out;
  const size_t nshell  = std::accumulate(geom_->atoms().begin(), geom_->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  out.reserve(nshell*nshell);

  // TODO perhaps we could reduce operation by a factor of 2
  int cnt = 0;
  int iatom0 = 0;
  auto oa0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = oa1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          // static distribution since this is cheap
          if (cnt++ % mpi__->size() != mpi__->rank()) continue;

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          vector<int> atom = {iatom0, iatom1};
          vector<int> offset_ = {*o0, *o1};

          out.push_back(make_shared<TaskType>(input, atom, offset_, nmat, kmat, omat, this));
        }
      }
    }
  }

  // if finite nucleus are used, we need to insert additional tasks that are omitted from the 1e part
  if (geom_->has_finite_nucleus()) {
    vector<shared_ptr<GradTask>> task0 = contract_grad1e_fnai(nmat);
    out.insert(out.end(), task0.begin(), task0.end());
  }

  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_graddkh1e(array<shared_ptr<const Matrix>, 4> den) {
  auto geom = make_shared<Molecule>(*geom_);
  geom = geom->uncontract();
  vector<shared_ptr<GradTask>> out;
  const size_t nshell  = std::accumulate(geom->atoms().begin(), geom->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  out.reserve(nshell*nshell);

  // TODO perhaps we could reduce operation by a factor of 2
  int cnt = 0;
  int iatom0 = 0;
  auto oa0 = geom->offsets().begin();
  for (auto a0 = geom->atoms().begin(); a0 != geom->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = geom->offsets().begin();
    for (auto a1 = geom->atoms().begin(); a1 != geom->atoms().end(); ++a1, ++oa1, ++iatom1) {

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = oa1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          // static distribution since this is cheap
          if (cnt++ % mpi__->size() != mpi__->rank()) continue;

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          vector<int> atom = {iatom0, iatom1};
          vector<int> offset_ = {*o0, *o1};

          out.push_back(make_shared<GradTask1d>(input, atom, offset_, den, this));
        }
      }
    }
  }

  return out;
}


// TODO make a generic code to merge with grad1e (variadic templete? vector?)
vector<shared_ptr<GradTask>> GradEval_base::contract_gradsmall1e(array<shared_ptr<const Matrix>,6> rmat) {
  vector<shared_ptr<GradTask>> out;
  const size_t nshell  = std::accumulate(geom_->atoms().begin(), geom_->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  out.reserve(nshell*nshell);

  // TODO perhaps we could reduce operation by a factor of 2
  int cnt = 0;
  int iatom0 = 0;
  auto oa0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = oa1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          // static distribution since this is cheap
          if (cnt++ % mpi__->size() != mpi__->rank()) continue;

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          vector<int> atom = {iatom0, iatom1};
          vector<int> offset_ = {*o0, *o1};

          out.push_back(make_shared<GradTask1r>(input, atom, offset_, rmat, this));
        }
      }
    }
  }
  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_grad2e(const array<shared_ptr<const DFDist>,6> o, const shared_ptr<const Geometry> geom) {
  shared_ptr<const Geometry> cgeom = geom ? geom : geom_;
  vector<shared_ptr<GradTask>> out;
  const size_t nshell  = std::accumulate(cgeom->atoms().begin(), cgeom->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  const size_t nshell2  = std::accumulate(cgeom->aux_atoms().begin(), cgeom->aux_atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });

  out.reserve(nshell*nshell*nshell2);

  int iatom0 = 0;
  auto oa0 = cgeom->offsets().begin();
  for (auto a0 = cgeom->atoms().begin(); a0 != cgeom->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = cgeom->offsets().begin();
    for (auto a1 = cgeom->atoms().begin(); a1 != cgeom->atoms().end(); ++a1, ++oa1, ++iatom1) {
      int iatom2 = 0;
      auto oa2 = cgeom->aux_offsets().begin();
      for (auto a2 = cgeom->aux_atoms().begin(); a2 != cgeom->aux_atoms().end(); ++a2, ++oa2, ++iatom2) {
        if ((*a2)->dummy()) continue;
        // dummy shell
        auto b3 = make_shared<const Shell>((*a2)->shells().front()->spherical());

        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = oa1->begin();
          for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {
            auto o2 = oa2->begin();
            for (auto b2 = (*a2)->shells().begin(); b2 != (*a2)->shells().end(); ++b2, ++o2) {
              tuple<size_t, size_t> info = o[0]->adist_now()->locate(*o2);
              if (get<0>(info) != mpi__->rank()) continue;

              array<shared_ptr<const Shell>,4> input = {{b3, *b2, *b1, *b0}};
              vector<int> atoms = {iatom0, iatom1, iatom2};
              vector<int> offs = {*o0, *o1, *o2};

              out.push_back(make_shared<GradTask3r>(input, atoms, offs, o, this));
            }
          }
        }

      }
    }
  }
  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_grad1e_fnai(const array<shared_ptr<const Matrix>,6> o, const shared_ptr<const Geometry> geom) {
  shared_ptr<const Geometry> cgeom = geom ? geom : geom_;
  vector<shared_ptr<GradTask>> out;
  const size_t nshell = std::accumulate(cgeom->atoms().begin(), cgeom->atoms().end(), 0, [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  const size_t nfatom = std::accumulate(cgeom->atoms().begin(), cgeom->atoms().end(), 0, [](const int& i, const shared_ptr<const Atom>& o) { return i+(o->finite_nucleus() ? 1 : 0); });
  out.reserve(nshell*nshell*nfatom);

  // loop over atoms
  int iatomf = -1;
  int cnt = 0;
  for (auto& af : cgeom->atoms()) {
    ++iatomf;
    if (!af->finite_nucleus()) continue;
    // construct nuclear shell
    const double fac = - af->atom_charge()*pow(af->atom_exponent()/pi__, 1.5);
    auto b2 = make_shared<Shell>(af->spherical(), af->position(), 0, vector<double>{af->atom_exponent()}, vector<vector<double>>{{fac}}, vector<pair<int,int>>{make_pair(0,1)});

    // dummy shell
    auto b3 = make_shared<Shell>(af->spherical());

    int iatom0 = 0;
    auto oa0 = cgeom->offsets().begin();
    for (auto a0 = cgeom->atoms().begin(); a0 != cgeom->atoms().end(); ++a0, ++oa0, ++iatom0) {
      int iatom1 = 0;
      auto oa1 = cgeom->offsets().begin();
      for (auto a1 = cgeom->atoms().begin(); a1 != cgeom->atoms().end(); ++a1, ++oa1, ++iatom1) {
        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = oa1->begin();
          for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

            // static distribution since this is cheap
            if (cnt++ % mpi__->size() != mpi__->rank()) continue;

            array<shared_ptr<const Shell>,4> input = {{b3, b2, *b1, *b0}};
            vector<int> atoms = {iatom0, iatom1, iatomf};
            vector<int> offs = {*o0, *o1, 0};

            out.push_back(make_shared<GradTask1rf>(input, atoms, offs, o, this));
          }

        }
      }
    }
  }
  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_grad2e(const shared_ptr<const DFDist> o, const shared_ptr<const Geometry> geom) {
  shared_ptr<const Geometry> cgeom = geom ? geom : geom_;
  vector<shared_ptr<GradTask>> out;
  const size_t nshell  = std::accumulate(cgeom->atoms().begin(), cgeom->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  const size_t nshell2  = std::accumulate(cgeom->aux_atoms().begin(), cgeom->aux_atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });

  out.reserve(nshell*(nshell+1)*nshell2/2);

  // loop over atoms (using symmetry b0 <-> b1)
  int iatom0 = 0;
  auto oa0 = cgeom->offsets().begin();
  for (auto a0 = cgeom->atoms().begin(); a0 != cgeom->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = iatom0;
    auto oa1 = oa0;
    for (auto a1 = a0; a1 != cgeom->atoms().end(); ++a1, ++oa1, ++iatom1) {
      int iatom2 = 0;
      auto oa2 = cgeom->aux_offsets().begin();
      for (auto a2 = cgeom->aux_atoms().begin(); a2 != cgeom->aux_atoms().end(); ++a2, ++oa2, ++iatom2) {
        if ((*a2)->dummy()) continue;
        // dummy shell
        auto b3 = make_shared<const Shell>((*a2)->shells().front()->spherical());

        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = a0!=a1 ? oa1->begin() : o0;
          for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {
            auto o2 = oa2->begin();
            for (auto b2 = (*a2)->shells().begin(); b2 != (*a2)->shells().end(); ++b2, ++o2) {
              tuple<size_t, size_t> info = o->adist_now()->locate(*o2);
              if (get<0>(info) != mpi__->rank()) continue;

              array<shared_ptr<const Shell>,4> input = {{b3, *b2, *b1, *b0}};
              vector<int> atoms = {iatom0, iatom1, iatom2};
              vector<int> offs = {*o0, *o1, *o2};

              out.push_back(make_shared<GradTask3>(input, atoms, offs, o, this));
            }
          }
        }

      }
    }
  }
  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_grad2e_2index(const shared_ptr<const Matrix> den, const shared_ptr<const Geometry> geom) {
  shared_ptr<const Geometry> cgeom = geom ? geom : geom_;
  vector<shared_ptr<GradTask>> out;
  const size_t nshell2  = std::accumulate(cgeom->aux_atoms().begin(), cgeom->aux_atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  out.reserve(nshell2*(nshell2+1)/2);

  // using symmetry (b0 <-> b1)
  int cnt = 0;
  int iatom0 = 0;
  auto oa0 = cgeom->aux_offsets().begin();
  for (auto a0 = cgeom->aux_atoms().begin(); a0 != cgeom->aux_atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = iatom0;
    auto oa1 = oa0;
    for (auto a1 = a0; a1 != cgeom->aux_atoms().end(); ++a1, ++oa1, ++iatom1) {
      if ((*a0)->dummy() || (*a1)->dummy()) continue;

      // dummy shell
      auto b3 = make_shared<const Shell>((*a0)->shells().front()->spherical());

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = a0!=a1 ? oa1->begin() : o0;
        for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          // static distribution since this is cheap
          if (cnt++ % mpi__->size() != mpi__->rank()) continue;

          array<shared_ptr<const Shell>,4> input = {{*b1, b3, *b0, b3}};
          vector<int> atoms = {iatom0, iatom1};
          vector<int> offs = {*o0, *o1};

          out.push_back(make_shared<GradTask2>(input, atoms, offs, den, this));
        }
      }
    }
  }
  return out;
}


vector<shared_ptr<GradTask>> GradEval_base::contract_grad1e_fnai(const shared_ptr<const Matrix> nmat) {
  vector<shared_ptr<GradTask>> out;
  const size_t nshell = std::accumulate(geom_->atoms().begin(), geom_->atoms().end(), 0, [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  const size_t nfatom = std::accumulate(geom_->atoms().begin(), geom_->atoms().end(), 0, [](const int& i, const shared_ptr<const Atom>& o) { return i+(o->finite_nucleus() ? 1 : 0); });
  out.reserve(nshell*nshell*nfatom);

  int cnt = 0;

  // loop over finite nucleus
  int iatomf = -1;
  for (auto& af : geom_->atoms()) {
    ++iatomf;
    // skip if this atom is not finite nucleus
    if (!af->finite_nucleus()) continue;

    // construct nuclear shell
    const double fac = - af->atom_charge()*pow(af->atom_exponent()/pi__, 1.5);
    auto nshell = make_shared<Shell>(af->spherical(), af->position(), 0, vector<double>{af->atom_exponent()}, vector<vector<double>>{{fac}}, vector<pair<int,int>>{make_pair(0,1)});

    // construct dummy shell
    auto dummy = make_shared<const Shell>(af->spherical());

    int iatom0 = 0;
    auto oa0 = geom_->offsets().begin();
    for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++oa0, ++iatom0) {
      int iatom1 = 0;
      auto oa1 = geom_->offsets().begin();
      for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++oa1, ++iatom1) {

        auto o0 = oa0->begin();
        for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
          auto o1 = oa1->begin();
          for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

            // static distribution since this is cheap
            if (cnt++ % mpi__->size() != mpi__->rank()) continue;

            array<shared_ptr<const Shell>,4> input = {{dummy, nshell, *b1, *b0}};
            vector<int> atom = {iatom0, iatom1, iatomf};
            vector<int> offset_ = {*o0, *o1, 0};

            out.push_back(make_shared<GradTask1f>(input, atom, offset_, nmat, this));
          }
        }
      }
    }
  }
  return out;
}

// explicit instantiation of the template functions
template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1>(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1s>(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1>(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k,
                                                                                 const std::shared_ptr<const Matrix> o);
template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1s>(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k,
                                                                                  const std::shared_ptr<const Matrix> o);
