//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks13.cc
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

#include <src/smith/caspt2/CASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task600::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I760
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task601::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I760
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, a2, c4);
}

void Task602::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I764
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task603::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I764
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, a2, c4);
}

void Task604::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
        // tensor label: I768
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task605::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I768
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, a2, c1, c4);
}

void Task606::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
        // tensor label: I772
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task607::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I772
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<1,0,2,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, c1, a4);
}

void Task608::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I776
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task609::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I776
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<2,0,1,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, x0, a2, c1, c4);
}

void Task610::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I780
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task611::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I780
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<1,0,2,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, c1, a4);
}

void Task612::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
    // tensor label: I800
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task613::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I800
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,2,1,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->add_block(odata, x0, a4);
}

void Task614::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I800
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->add_block(odata, x0, a4);
}

void Task615::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size());
    // tensor label: I943
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task616::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I943
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,3,2,1,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), a2.size(), a4.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*c3.size()*c1.size(), i1data_sorted, a4.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->add_block(odata, a2, x1);
}

void Task617::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size());
    // tensor label: I947
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task618::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I947
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->add_block(odata, a4, x1);
}

void Task619::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x0);
        dscal_(x0.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task620::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x0);
        dscal_(x0.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task621::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task622::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x0, x1);
}

void Task623::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0, c1, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size(), c1.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task624::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I759
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x0, x1);
}

void Task625::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I698
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma270
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x1, x0, x2)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
          // tensor label: I783
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x3, x2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x3, x2)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x3.size(), x2.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x3.size()*x2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*x3.size()*x2.size(), i1data_sorted, x0.size()*x1.size()*x3.size()*x2.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task626::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I783
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I784
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x1, x3, x2);
}

void Task627::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I784
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, a2, x1);
}

void Task628::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I783
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x0, x1, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x0, x1, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), x0.size(), x1.size(), a1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x0, x1, x3, x2);
}

void Task629::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I698
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma276
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x3, x2, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x1, x3, x2, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
              // tensor label: I807
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, x5, x4, x3)]);
              sort_indices<3,4,0,5,2,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
              dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
                     1.0, i0data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(), i1data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task630::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I807
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    // tensor label: I808
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, x0, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), x0.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task631::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I808
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c2, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, c2, x0, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2, x0, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2, x0, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
    sort_indices<1,0,2,3,0,1,-2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*c2.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), x0.size(), c2.size(), x1.size());
  out()->add_block(odata, x1, c2, x0, x2);
}

void Task632::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I698
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: Gamma277
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x3, x2, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x1, x3, x2, x0)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
          // tensor label: I811
          std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
          dgemm_("T", "N", ci0.size(), 1, x2.size()*x3.size()*x1.size()*x0.size(),
                 1.0, i0data_sorted, x2.size()*x3.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x3.size()*x1.size()*x0.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task633::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I811
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I812
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, c2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, c2, x3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), c2.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), a1.size()*c2.size(),
             1.0, i0data_sorted, a1.size()*c2.size(), i1data_sorted, a1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->add_block(odata, x2, x3, x1, x0);
}

void Task634::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I812
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: f1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", a1.size()*c2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a1, c2, x3);
}

void Task635::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I811
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3, x2, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x2, x3, x1, x0);
}

void Task636::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I698
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: Gamma279
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x1, x4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x1, x4)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
          // tensor label: I819
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x5, x4)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x5.size(), x4.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x5.size()*x4.size(), i1data_sorted, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task637::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I819
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, x5, x4);
}

void Task638::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I698
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma280
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x0, x1, x2)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
          // tensor label: I822
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x3, x2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x3, x2)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x3.size(), x2.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x3.size()*x2.size(), i1data_sorted, x1.size()*x0.size()*x3.size()*x2.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task639::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), x2.size());
      // tensor label: I823
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, x0, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, x0, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), x0.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a1.size()*c3.size(),
             1.0, i0data_sorted, a1.size()*c3.size(), i1data_sorted, a1.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task640::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  // tensor label: I823
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a1, x0, c3)]);
  std::fill_n(odata.get(), out()->get_size(x1, a1, x0, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a1, x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a1, x0, c3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
    sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
    dgemm_("T", "N", c3.size(), x1.size()*a1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size(), a1.size(), x1.size());
  out()->add_block(odata, x1, a1, x0, c3);
}

void Task641::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I827
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, x0, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, x0, a3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), x0.size(), a3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task642::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  // tensor label: I827
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c2, x0, a3)]);
  std::fill_n(odata.get(), out()->get_size(x1, c2, x0, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2, x0, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2, x0, a3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
    sort_indices<1,0,2,3,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
    dgemm_("T", "N", a3.size(), x1.size()*c2.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), x1.size());
  out()->add_block(odata, x1, c2, x0, a3);
}

void Task643::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I858
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a1, c2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), c2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), a1.size()*c2.size(),
             1.0, i0data_sorted, a1.size()*c2.size(), i1data_sorted, a1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task644::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I858
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: f1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a1.size()*c2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*c2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), c2.size(), x2.size());
  out()->add_block(odata, x2, x3, a1, c2);
}

void Task645::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I985
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task646::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I985
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    dscal_(x1.size()*c2.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,-2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->add_block(odata, c2, a1, x0, x1);
}

void Task647::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I985
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x1.size(), c2.size()*a1.size()*x0.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, a1, x0, x1);
}

void Task648::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task649::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

#endif
