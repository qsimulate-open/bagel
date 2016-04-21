//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks12.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task550::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I645
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());
    // tensor label: I646
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    dgemm_("T", "N", c2.size(), x1.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size(), a3.size(), a1.size());
  out()->add_block(odata, x1, a3, a1, c2);
}

void Task551::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);
  // tensor label: I646
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, c4, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, c4, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a3.size(), c4.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c4, a3);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x1.size(), a1.size(), c4.size(), a3.size());
  }
  out()->add_block(odata, x1, a3, c4, a1);
}

void Task552::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I645
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I654
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    dgemm_("T", "N", a1.size(), x1.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, a1.size(), x1.size(), c2.size(), a3.size());
  out()->add_block(odata, x1, a3, a1, c2);
}

void Task553::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I654
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4, c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4, c2, a3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), a4.size(), c2.size(), a3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a3, c2, a4);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x1.size(), a3.size(), c2.size(), a4.size());
  }
  out()->add_block(odata, x1, a4, c2, a3);
}

void Task554::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I645
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a3, a1, c2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I662
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    dgemm_("T", "N", a3.size(), x1.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size(), c2.size(), a1.size());
  out()->add_block(odata, x1, a3, a1, c2);
}

void Task555::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I662
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a4.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, a4);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a4.size());
  }
  out()->add_block(odata, x1, a4, c2, a1);
}

void Task556::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), x0.size());
      // tensor label: I834
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task557::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I834
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I835
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      dgemm_("T", "N", 1, c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size());
  out()->add_block(odata, c2, a1);
}

void Task558::Task_local::compute() {
  const Index c2 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a1 = b(3);
  // tensor label: I835
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a4, c3, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a4, c3, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c2.size(), a4.size(), c3.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a1, c3, a4);
    sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c2.size(), a1.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c2, a4, c3, a1);
}

void Task559::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I888
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, c1.size()*a2.size(),
             1.0, i0data_sorted, c1.size()*a2.size(), i1data_sorted, c1.size()*a2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task560::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I888
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I889
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, c1, a2);
}

void Task561::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I889
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->add_block(odata, c1, a4, c3, a2);
}

void Task562::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I939
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task563::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I939
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,1,0,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->add_block(odata, a4, c1);
}

void Task564::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I943
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task565::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I943
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, a2, c1);
}

void Task566::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x0.size(), x1.size());
      // tensor label: I947
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task567::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I947
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,1,0,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->add_block(odata, a4, c1);
}

void Task568::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I951
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task569::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I951
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
      // tensor label: l2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, a2, c1);
}

void Task570::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I961
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task571::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I961
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, a2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), a2.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<3,0,1,2,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x0.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task572::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I961
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), a4.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<1,0,3,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x0.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task573::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), c4.size());
    // tensor label: I993
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), x1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task574::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I993
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c4)]);
  std::fill_n(odata.get(), out()->get_size(x1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<2,3,1,0,0,1,-4,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", c4.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size());
  out()->add_block(odata, x1, c4);
}

void Task575::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I993
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c4)]);
  std::fill_n(odata.get(), out()->get_size(x1, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<2,1,3,0,0,1,2,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", c4.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x1.size());
  out()->add_block(odata, x1, c4);
}

void Task576::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: l2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I1007
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, a1, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, a1, c2)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), a1.size(), c2.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*a1.size()*c2.size(),
               1.0, i0data_sorted, a3.size()*a1.size()*c2.size(), i1data_sorted, a3.size()*a1.size()*c2.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task577::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1007
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());
    // tensor label: I1008
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());
    dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
  out()->add_block(odata, x0, a3, a1, c2);
}

void Task578::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);
  // tensor label: I1008
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata.get(), out()->get_size(x0, a3, c4, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c4, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c4.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a1, c4, a3);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x0.size(), a1.size(), c4.size(), a3.size());
  }
  out()->add_block(odata, x0, a3, c4, a1);
}

void Task579::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1007
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I1016
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());
    dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
  out()->add_block(odata, x0, a3, a1, c2);
}

void Task580::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I1016
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4, c2, a3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x0.size(), a4.size(), c2.size(), a3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a3, c2, a4);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x0.size(), a3.size(), c2.size(), a4.size());
  }
  out()->add_block(odata, x0, a4, c2, a3);
}

void Task581::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1007
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, a1, c2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I1024
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());
    dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, x0, a3, a1, c2);
}

void Task582::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1024
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x0.size(), a4.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a1, c2, a4);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x0.size(), a1.size(), c2.size(), a4.size());
  }
  out()->add_block(odata, x0, a4, c2, a1);
}

void Task583::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
        dscal_(a3.size()*c2.size()*a1.size()*x0.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task584::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
        dscal_(a3.size()*c2.size()*a1.size()*x0.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,-2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task585::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        dscal_(a3.size()*c2.size()*a1.size()*x1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task586::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I471
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, c2, a3);
        dscal_(a3.size()*c2.size()*a1.size()*x1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,-2,1>(i1data, i1data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task587::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma161
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I521
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, x5, x4, x3)]);
              sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
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

void Task588::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I521
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: I522
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task589::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I522
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, a2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, a2, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x0, x1)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*x0.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), a2.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, a2, x2);
}

void Task590::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I521
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I908
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x1, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", x3.size()*x4.size()*x5.size(), x0.size()*x1.size()*x2.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size()*x4.size()*x5.size());
  }
  sort_indices<4,3,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task591::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I908
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x0, x1)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
    dgemm_("T", "N", x2.size(), a1.size()*x0.size()*x1.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), a1.size(), x0.size(), x1.size());
  out()->add_block(odata, a1, x0, x1, x2);
}

void Task592::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma166
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x4, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x3, x4, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
              // tensor label: I541
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x2, x1, x0)]);
              sort_indices<0,5,2,1,3,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", ci0.size(), 1, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
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

void Task593::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I541
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I542
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task594::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I542
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, a1, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x5, a1, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, a1, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, a1, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    dgemm_("T", "N", x3.size(), x5.size()*a1.size()*x4.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x5.size(), a1.size(), x4.size());
  out()->add_block(odata, x5, a1, x4, x3);
}

void Task595::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I541
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I830
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, x5, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, x5, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x5.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<4,3,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task596::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I830
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, a1, x5, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, a1, x5, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, a1, x5, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x5, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    dgemm_("T", "N", x3.size(), x4.size()*a1.size()*x5.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x5.size(), a1.size(), x4.size());
  out()->add_block(odata, x4, a1, x5, x3);
}

void Task597::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x0 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma167
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
              // tensor label: I545
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x2, x1, x0)]);
              sort_indices<0,1,2,5,3,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", ci0.size(), 1, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x2.size()*x1.size()*x0.size(),
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

void Task598::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I545
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: l2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I546
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task599::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I546
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x5, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x5, x4)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
    dgemm_("T", "N", x3.size(), a1.size()*x5.size()*x4.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x5.size(), x4.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

#endif
