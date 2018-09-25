//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector_impl.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ci/ras/civector_impl.h>

using namespace std;
using namespace bagel;


template<typename DataType>
void RASCivector_impl<DataType>::spin_decontaminate(const double thresh) {
  const int nspin = det_->nspin();
  const int max_spin = det_->nelea() + det_->neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  auto S2 = spin_();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while(fabs(actual_expectation - pure_expectation) > thresh) {
    if (k > max_spin) { this->print(0.05); throw std::runtime_error("Spin decontamination failed."); }

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin_();
    actual_expectation = dot_product(*S2);

    k += 2;
  }
}


template<typename DataType>
double RASCivector_impl<DataType>::normalize() {
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  scale(DataType(scal));
  return norm;
}


template<typename DataType>
void RASCivector_impl<DataType>::print(const double thr) const {
  // multimap sorts elements so that they will be shown in the descending order in magnitude
  multimap<double, tuple<DataType, bitset<nbit__>, bitset<nbit__>>> tmp;
  for (auto& iblock : blocks_) {
    if (!iblock) continue;
    double* i = iblock->data();
    for (auto& ia : *iblock->stringsa()) {
      for (auto& ib : *iblock->stringsb()) {
        if (abs(*i) > thr)
          tmp.emplace(-abs(*i), make_tuple(*i, ia, ib));
        ++i;
      }
    }
  }
  for (auto& i : tmp)
    cout << "       " << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0))
              << "-" << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0), det_->ras(0)+det_->ras(1))
              << "-" << print_bit(get<1>(i.second), get<2>(i.second), det_->ras(0)+det_->ras(1), det_->norb())
              << "  " << setprecision(10) << setw(15) << get<0>(i.second) << endl;
}


template<typename DataType>
void RASCivector_impl<DataType>::synchronize(const int root) {
#ifdef HAVE_MPI_H
  mpi__->broadcast(data(), size(), root);
#endif /* HAVE_MPI_H */
}


template class bagel::RASCivector_impl<double>;

