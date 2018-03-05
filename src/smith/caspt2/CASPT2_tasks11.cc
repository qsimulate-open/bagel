//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks11.cc
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

#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task500::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I660
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task501::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I662
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    sort_indices<0,3,2,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x0.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task502::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I662
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I663
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, c3, a2, c1, x0);
}

void Task503::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I663
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->add_block(odata, c3, a2, c1, x1);
}

void Task504::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  {
    // tensor label: I666
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, c2, x0, x1);
    sort_indices<1,3,2,0,1,1,1,1>(i0data, odata, a1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task505::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I666
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task506::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I666
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task507::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, c1, a2), 0.0);
  {
    // tensor label: I670
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1, x1, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x0, x1, c1, a2);
}

void Task508::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I670
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I671
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, a2, c1, x1, x0);
}

void Task509::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I671
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c1, x2);
}

void Task510::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I674
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task511::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I674
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma59
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task512::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, a2), 0.0);
  {
    // tensor label: I676
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->add_block(odata, c3, a4, c1, a2);
}

void Task513::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I676
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

void Task514::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3, x0, a1), 0.0);
  {
    // tensor label: I678
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, a1, x0);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a3.size(), c2.size(), a1.size(), x0.size());
  }
  out()->add_block(odata, c2, a3, x0, a1);
}

void Task515::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I678
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I679
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->add_block(odata, a3, c2, a1, x0);
}

void Task516::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I679
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x1, a3, c2, a1);
}

void Task517::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, x0, a1), 0.0);
  {
    // tensor label: I682
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a2, x0, x1);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, a1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x1, a2, x0, a1);
}

void Task518::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I682
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a2)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->add_block(odata, a1, a2, x0, x1);
}

void Task521::Task_local::compute() {
  const Index x0 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: I685
  std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x0)]);
  sort_indices<3,0,2,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
  // tensor label (calculated on-the-fly): Gamma248
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
  if (x1 == x5 && x0 == x4) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            o1data[ix3+x3.size()*(ix2)] +=
              -2.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix4+x0.size()*(ix5+x5.size()*(ix5+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x2);
  }
  if (x1 == x4 && x0 == x5) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            o1data[ix3+x3.size()*(ix2)] +=
              4.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix5+x0.size()*(ix5+x5.size()*(ix4+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x2);
  }
  // rdm0 merged ci derivative case
  if (x0 == x2 && x3 == x5 && x1 == x4) {
    std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
    std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          o0data[0] +=
            4.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix4+x1.size()*(ix4)))];
        }
      }
    }
    out(0)->add_block(o0data);
  }
  if (x1 == x4 && x0 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x5)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix3+x3.size()*(ix5)] +=
              -2.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix4+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x5);
  }
  // rdm0 merged ci derivative case
  if (x0 == x2 && x3 == x4 && x1 == x5) {
    std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
    std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          o0data[0] +=
            -2.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix5+x1.size()*(ix4)))];
        }
      }
    }
    out(0)->add_block(o0data);
  }
  if (x1 == x5 && x0 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x4)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix3+x3.size()*(ix4)] +=
              1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix5+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x4);
  }
  if (x3 == x4 && x0 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x5)]);
    std::fill_n(o1data.get(), out(1)->get_size(x1, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix1+x1.size()*(ix5)] +=
              1.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x1, x5);
  }
  if (x3 == x5 && x0 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x4)]);
    std::fill_n(o1data.get(), out(1)->get_size(x1, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix1+x1.size()*(ix4)] +=
              -2.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x1, x4);
  }
  if (x0 == x2) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x3, x5, x1, x4)]);
    std::fill_n(o2data.get(), out(2)->get_size(x3, x5, x1, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              o2data[ix3+x3.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))] +=
                1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix2+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x3, x5, x1, x4);
  }
  // rdm0 merged ci derivative case
  if (x0 == x4 && x1 == x2 && x3 == x5) {
    std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
    std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          o0data[0] +=
            -2.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix4+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
        }
      }
    }
    out(0)->add_block(o0data);
  }
  if (x0 == x4 && x1 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x5)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix3+x3.size()*(ix5)] +=
              1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix4+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x5);
  }
  if (x3 == x5 && x0 == x4) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x1, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            o1data[ix1+x1.size()*(ix2)] +=
              1.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix4+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x1, x2);
  }
  if (x0 == x4) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x1, x5, x3, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x1, x5, x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              o2data[ix1+x1.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2)))] +=
                1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix4+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x1, x5, x3, x2);
  }
  // rdm0 merged ci derivative case
  if (x3 == x4 && x0 == x5 && x1 == x2) {
    std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
    std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          o0data[0] +=
            4.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix5+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
        }
      }
    }
    out(0)->add_block(o0data);
  }
  if (x0 == x5 && x1 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x4)]);
    std::fill_n(o1data.get(), out(1)->get_size(x3, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix3+x3.size()*(ix4)] +=
              -2.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix5+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x3, x4);
  }
  if (x3 == x4 && x0 == x5) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x1, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            o1data[ix1+x1.size()*(ix2)] +=
              -2.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix5+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x1, x2);
  }
  if (x0 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x1, x4, x3, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x1, x4, x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o2data[ix1+x1.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2)))] +=
                -2.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix5+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x1, x4, x3, x2);
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x5)]);
    std::fill_n(o1data.get(), out(1)->get_size(x0, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix0+x0.size()*(ix5)] +=
              -2.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x0, x5);
  }
  if (x3 == x5 && x1 == x2) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x4)]);
    std::fill_n(o1data.get(), out(1)->get_size(x0, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            o1data[ix0+x0.size()*(ix4)] +=
              1.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x0, x4);
  }
  if (x1 == x2) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x0, x5, x3, x4)]);
    std::fill_n(o2data.get(), out(2)->get_size(x0, x5, x3, x4), 0.0);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              o2data[ix0+x0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4)))] +=
                1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix2+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x0, x5, x3, x4);
  }
  if (x3 == x5 && x1 == x4) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x0, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            o1data[ix0+x0.size()*(ix2)] +=
              -2.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix4+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x0, x2);
  }
  if (x1 == x4) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x0, x5, x3, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x0, x5, x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              o2data[ix0+x0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2)))] +=
                -2.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix4+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x0, x5, x3, x2);
  }
  if (x3 == x4 && x1 == x5) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x2)]);
    std::fill_n(o1data.get(), out(1)->get_size(x0, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            o1data[ix0+x0.size()*(ix2)] +=
              1.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix5+x1.size()*(ix4)))];
          }
        }
      }
    }
    out(1)->add_block(o1data, x0, x2);
  }
  if (x1 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x0, x4, x3, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x0, x4, x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o2data[ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2)))] +=
                1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix5+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x0, x4, x3, x2);
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x0, x5, x1, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x0, x5, x1, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              o2data[ix0+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix2)))] +=
                1.0 * fdata[ix4+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x0, x5, x1, x2);
  }
  if (x3 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x1, x4, x0, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x1, x4, x0, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o2data[ix1+x1.size()*(ix4+x4.size()*(ix0+x0.size()*(ix2)))] +=
                1.0 * fdata[ix5+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x1, x4, x0, x2);
  }
  {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x0, x5, x1, x4, x3, x2)]);
    std::fill_n(o3data.get(), out(3)->get_size(x0, x5, x1, x4, x3, x2), 0.0);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
                o3data[ix0+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2)))))] +=
                  1.0 * fdata[ix3+x3.size()*(ix2)] * i0data_sorted[ix0+x0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix4)))];
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x0, x5, x1, x4, x3, x2);
  }
}

void Task522::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I685
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x5, c2, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
      sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      dgemm_("T", "N", x1.size()*x0.size(), x5.size()*x4.size(), c1.size()*c2.size(),
             1.0, i0data_sorted, c1.size()*c2.size(), i1data_sorted, c1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x5.size(), x4.size());
  out()->add_block(odata, x5, x4, x1, x0);
}

void Task523::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, x3, x2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0, x3, x2)]);
  sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size(), x3.size(), x2.size());
  // tensor label (calculated on-the-fly): Gamma249
  {
    // rdm0 non-merged ci derivative case
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
      std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          o0data[0] +=
            -2.0 * i0data_sorted[ix2+x0.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2)))];
        }
      }
      out(0)->add_block(o0data);
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            o1data[ix1+x1.size()*(ix3)] +=
              1.0 * i0data_sorted[ix2+x0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x1, x3);
    }
  }
  {
    // rdm0 non-merged ci derivative case
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
      std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          o0data[0] +=
            4.0 * i0data_sorted[ix3+x0.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2)))];
        }
      }
      out(0)->add_block(o0data);
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            o1data[ix1+x1.size()*(ix2)] +=
              -2.0 * i0data_sorted[ix3+x0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x1, x2);
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x0, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            o1data[ix0+x0.size()*(ix3)] +=
              -2.0 * i0data_sorted[ix0+x0.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x0, x3);
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x0, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            o1data[ix0+x0.size()*(ix2)] +=
              1.0 * i0data_sorted[ix0+x0.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x0, x2);
    }
  }
  {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x0, x3, x1, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x0, x3, x1, x2), 0.0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data_sorted, o2data, x0.size(), x3.size(), x1.size(), x2.size());
    out(2)->add_block(o2data, x0, x3, x1, x2);
  }
}

void Task524::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      // tensor label: I689
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c1.size()*c3.size(),
             1.0, i0data_sorted, c1.size()*c3.size(), i1data_sorted, c1.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task525::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I689
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, c1, c3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, c3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, c2, x1)]);
    sort_indices<2,0,1,3,0,1,-8,1>(i1data, i1data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
    dgemm_("T", "N", c3.size(), x1.size()*x0.size()*c1.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, c1, c3);
}

void Task526::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x3, c2, x2);
      dscal_(c1.size()*x3.size()*c2.size()*x2.size(), e0_, i1data.get(), 1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
      sort_indices<2,0,1,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x3.size()*x2.size(), c1.size()*c2.size(),
             1.0, i0data_sorted, c1.size()*c2.size(), i1data_sorted, c1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task527::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, c2, x1)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), c2.size()*c1.size(),
             1.0, i0data_sorted, c2.size()*c1.size(), i1data_sorted, c2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task528::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x3, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c2, x2)]);
      sort_indices<0,2,1,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), c2.size()*c1.size(),
             1.0, i0data_sorted, c2.size()*c1.size(), i1data_sorted, c2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x1, x0, x3, x2);
}

void Task529::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I692
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, x2, x5, x4, x3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0, x2, x5, x4, x3)]);
  sort_indices<3,4,1,5,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
  // tensor label (calculated on-the-fly): Gamma250
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x4)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x4), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix4)] +=
                -2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix2+x0.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x4);
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix3)] +=
                1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix2+x0.size()*(ix3+x3.size()*(ix4+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x3);
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x1, x3)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x1, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3)))] +=
                  1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix2+x0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x1, x3);
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x4)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x4), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix4)] +=
                4.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x4);
    }
  }
  {
    if (x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix2)] +=
                -2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix4+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x2);
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x1, x2)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x1, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix2)))] +=
                  -2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x1, x2);
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix3)] +=
                -2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x3);
    }
  }
  {
    if (x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix2)] +=
                1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x2);
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x3, x1, x2)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x3, x1, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))] +=
                  1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix1+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x3, x1, x2);
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x0, x3)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x0, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3)))] +=
                  -2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x0, x3);
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x0, x2)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x0, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix2)))] +=
                  1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x0, x2);
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x2, x0, x3)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x2, x0, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix2+x2.size()*(ix0+x0.size()*(ix3)))] +=
                  1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4+x1.size()*(ix2)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x2, x0, x3);
    }
  }
  {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x5, x4, x0, x3, x1, x2)]);
    std::fill_n(o3data.get(), out(3)->get_size(x5, x4, x0, x3, x1, x2), 0.0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data_sorted, o3data, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
    out(3)->add_block(o3data, x5, x4, x0, x3, x1, x2);
  }
}

void Task530::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I692
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
    // tensor label: I693
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task531::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I693
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, c2, x1)]);
    sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*x0.size()*c1.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, c1, x2);
}

void Task532::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I696
  std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
  sort_indices<2,0,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
  // tensor label (calculated on-the-fly): Gamma251
  {
    // rdm0 non-merged ci derivative case
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
      std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          o0data[0] +=
            -4.0 * i0data_sorted[ix3+x1.size()*(ix3+x3.size()*(ix2+x0.size()*(ix2)))];
        }
      }
      out(0)->add_block(o0data);
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            o1data[ix1+x1.size()*(ix3)] +=
              2.0 * i0data_sorted[ix1+x1.size()*(ix3+x3.size()*(ix2+x0.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x1, x3);
    }
  }
  {
    // rdm0 non-merged ci derivative case
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> o0data(new double[out(0)->get_size()]);
      std::fill_n(o0data.get(), out(0)->get_size(), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          o0data[0] +=
            2.0 * i0data_sorted[ix2+x1.size()*(ix3+x3.size()*(ix3+x0.size()*(ix2)))];
        }
      }
      out(0)->add_block(o0data);
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            o1data[ix1+x1.size()*(ix2)] +=
              -1.0 * i0data_sorted[ix1+x1.size()*(ix3+x3.size()*(ix3+x0.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x1, x2);
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x3)]);
      std::fill_n(o1data.get(), out(1)->get_size(x0, x3), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            o1data[ix0+x0.size()*(ix3)] +=
              -1.0 * i0data_sorted[ix2+x1.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x0, x3);
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x0, x2)]);
      std::fill_n(o1data.get(), out(1)->get_size(x0, x2), 0.0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            o1data[ix0+x0.size()*(ix2)] +=
              2.0 * i0data_sorted[ix3+x1.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2)))];
          }
        }
      }
      out(1)->add_block(o1data, x0, x2);
    }
  }
  {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x1, x3, x0, x2)]);
    std::fill_n(o2data.get(), out(2)->get_size(x1, x3, x0, x2), 0.0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data_sorted, o2data, x1.size(), x3.size(), x0.size(), x2.size());
    out(2)->add_block(o2data, x1, x3, x0, x2);
  }
}

void Task533::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I696
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: I697
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c2, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x3.size()*x2.size(), c1.size()*c2.size(),
             1.0, i0data_sorted, c1.size()*c2.size(), i1data_sorted, c1.size()*c2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task534::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I697
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a3, c2, x3)]);
    sort_indices<1,0,2,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    dgemm_("T", "N", x2.size(), c1.size()*c2.size()*x3.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), c2.size(), x3.size());
  out()->add_block(odata, c1, c2, x3, x2);
}

void Task535::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I696
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      // tensor label: I728
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c3.size()*c1.size(),
             1.0, i0data_sorted, c3.size()*c1.size(), i1data_sorted, c3.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task536::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I728
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, c1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<1,0,2,3,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, c1, x1);
}

void Task537::Task_local::compute() {
  const Index x2 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I700
  std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
  sort_indices<3,0,2,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  // tensor label (calculated on-the-fly): Gamma252
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x3, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix3+x3.size()*(ix0)] +=
                -2.0 * i0data_sorted[ix5+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix4+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x3, x0);
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x3, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x3, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix3+x3.size()*(ix0)] +=
                1.0 * i0data_sorted[ix4+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix5+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x3, x0);
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x2, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x2, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              o1data[ix2+x2.size()*(ix0)] +=
                1.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix5+x3.size()*(ix4+x4.size()*(ix4+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x2, x0);
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x2, x5, x3, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x2, x5, x3, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                o2data[ix2+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix0)))] +=
                  1.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix4+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x2, x5, x3, x0);
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x2, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x2, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              o1data[ix2+x2.size()*(ix0)] +=
                -2.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix4+x3.size()*(ix4+x4.size()*(ix5+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x2, x0);
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x3, x4, x2, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x3, x4, x2, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                o2data[ix3+x3.size()*(ix4+x4.size()*(ix2+x2.size()*(ix0)))] +=
                  1.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix5+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x3, x4, x2, x0);
    }
  }
  {
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix1+x1.size()*(ix0)] +=
                -2.0 * i0data_sorted[ix4+x2.size()*(ix5+x5.size()*(ix5+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x1, x0);
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x3, x5, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x3, x5, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix3+x3.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))] +=
                  1.0 * i0data_sorted[ix4+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x3, x5, x1, x0);
    }
  }
  {
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x1, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix1+x1.size()*(ix0)] +=
                4.0 * i0data_sorted[ix5+x2.size()*(ix5+x5.size()*(ix4+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x1, x0);
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x3, x4, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x3, x4, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix3+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))] +=
                  -2.0 * i0data_sorted[ix5+x2.size()*(ix5+x5.size()*(ix3+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x3, x4, x1, x0);
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x2, x5, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x2, x5, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                o2data[ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))] +=
                  -2.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix4+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x2, x5, x1, x0);
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x2, x4, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x2, x4, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                o2data[ix2+x2.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))] +=
                  1.0 * i0data_sorted[ix2+x2.size()*(ix5+x5.size()*(ix5+x3.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x2, x4, x1, x0);
    }
  }
  {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x2, x5, x3, x4, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x2, x5, x3, x4, x1, x0), 0.0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data_sorted, o3data, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
    out(3)->add_block(o3data, x2, x5, x3, x4, x1, x0);
  }
}

void Task538::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I700
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I701
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task539::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I701
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, c2, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, c2, x4)]);
    sort_indices<2,0,1,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
    dgemm_("T", "N", x3.size(), c1.size()*x5.size()*x4.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x5.size(), x4.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task540::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  const Index x4 = b(6);
  const Index x3 = b(7);
  // tensor label: I704
  std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x2, x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x2, x1, x0)]);
  sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
  // tensor label (calculated on-the-fly): Gamma253
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x4, x3);
  if (x2 == x6 && x1 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x0, x4, x3)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x0, x4, x3), 0.0);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix7+x7.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3)))] +=
                  -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix6+x2.size()*(ix5+x5.size()*(ix5+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x0, x4, x3);
  }
  if (x2 == x5 && x1 == x6) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x0, x4, x3)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x0, x4, x3), 0.0);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                o2data[ix7+x7.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3)))] +=
                  2.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix5+x2.size()*(ix5+x5.size()*(ix6+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x0, x4, x3);
  }
  if (x1 == x3 && x2 == x5 && x4 == x6) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x7, x0)]);
    std::fill_n(o1data.get(), out(1)->get_size(x7, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              o1data[ix7+x7.size()*(ix0)] +=
                2.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix5+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
            }
          }
        }
      }
    }
    out(1)->add_block(o1data, x7, x0);
  }
  if (x2 == x5 && x1 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x6, x4, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x6, x4, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix6+x6.size()*(ix4+x4.size()*(ix0)))] +=
                  2.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix5+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x6, x4, x0);
  }
  if (x4 == x5 && x2 == x6 && x1 == x3) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x7, x0)]);
    std::fill_n(o1data.get(), out(1)->get_size(x7, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              o1data[ix7+x7.size()*(ix0)] +=
                -1.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix6+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
            }
          }
        }
      }
    }
    out(1)->add_block(o1data, x7, x0);
  }
  if (x2 == x6 && x1 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x5, x4, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x5, x4, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0)))] +=
                  -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix6+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x5, x4, x0);
  }
  if (x4 == x5 && x1 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x6, x2, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x6, x2, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix0)))] +=
                  -1.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x6, x2, x0);
  }
  if (x4 == x6 && x1 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x0, x2, x5)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x0, x2, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix0+x0.size()*(ix2+x2.size()*(ix5)))] +=
                  -1.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x0, x2, x5);
  }
  if (x1 == x3) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x6, x2, x5, x4, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x6, x2, x5, x4, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                  o3data[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0)))))] +=
                    -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix3+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x6, x2, x5, x4, x0);
  }
  if (x4 == x6 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x7, x0)]);
    std::fill_n(o1data.get(), out(1)->get_size(x7, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
        for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix7+x7.size()*(ix0)] +=
                -1.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix5+x1.size()*(ix0)))))];
            }
          }
        }
      }
    }
    out(1)->add_block(o1data, x7, x0);
  }
  if (x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x6, x4, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x6, x4, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix7+x7.size()*(ix6+x6.size()*(ix4+x4.size()*(ix0)))] +=
                  -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix5+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x6, x4, x0);
  }
  if (x4 == x6 && x1 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x3, x2, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x3, x2, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix7+x7.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0)))] +=
                  -1.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x3, x2, x0);
  }
  if (x1 == x5) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x6, x4, x3, x2, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x6, x4, x3, x2, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  o3data[ix7+x7.size()*(ix6+x6.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0)))))] +=
                    -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x6, x4, x3, x2, x0);
  }
  if (x4 == x5 && x2 == x3 && x1 == x6) {
    std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x7, x0)]);
    std::fill_n(o1data.get(), out(1)->get_size(x7, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              o1data[ix7+x7.size()*(ix0)] +=
                2.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix6+x1.size()*(ix0)))))];
            }
          }
        }
      }
    }
    out(1)->add_block(o1data, x7, x0);
  }
  if (x2 == x3 && x1 == x6) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x0, x4, x5)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x0, x4, x5), 0.0);
    for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                o2data[ix7+x7.size()*(ix0+x0.size()*(ix4+x4.size()*(ix5)))] +=
                  -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix6+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x0, x4, x5);
  }
  if (x4 == x5 && x1 == x6) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x0, x2, x3)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x0, x2, x3), 0.0);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                o2data[ix7+x7.size()*(ix0+x0.size()*(ix2+x2.size()*(ix3)))] +=
                  -1.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix6+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x0, x2, x3);
  }
  if (x1 == x6) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x0, x2, x5, x4, x3)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x0, x2, x5, x4, x3), 0.0);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                  o3data[ix7+x7.size()*(ix0+x0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3)))))] +=
                    -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix6+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x0, x2, x5, x4, x3);
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x6, x1, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x6, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix6+x6.size()*(ix1+x1.size()*(ix0)))] +=
                  2.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x6, x1, x0);
  }
  if (x4 == x6 && x2 == x3) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x5, x1, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x5, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                o2data[ix7+x7.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))] +=
                  -1.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x5, x1, x0);
  }
  if (x2 == x3) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x6, x4, x5, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x6, x4, x5, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                  o3data[ix7+x7.size()*(ix6+x6.size()*(ix4+x4.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))] +=
                    -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix3+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x6, x4, x5, x1, x0);
  }
  if (x4 == x6 && x2 == x5) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x3, x1, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x3, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix7+x7.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))] +=
                  2.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix5+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x3, x1, x0);
  }
  if (x2 == x5) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x6, x4, x3, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x6, x4, x3, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  o3data[ix7+x7.size()*(ix6+x6.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))))] +=
                    2.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix5+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x6, x4, x3, x1, x0);
  }
  if (x4 == x5 && x2 == x6) {
    std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x7, x3, x1, x0)]);
    std::fill_n(o2data.get(), out(2)->get_size(x7, x3, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                o2data[ix7+x7.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))] +=
                  -1.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix6+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(2)->add_block(o2data, x7, x3, x1, x0);
  }
  if (x2 == x6) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x5, x4, x3, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x5, x4, x3, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                  o3data[ix7+x7.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))))] +=
                    -1.0 * fdata[ix4+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix6+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x5, x4, x3, x1, x0);
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x6, x2, x3, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x6, x2, x3, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  o3data[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))))] +=
                    -1.0 * fdata[ix5+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x6, x2, x3, x1, x0);
  }
  if (x4 == x6) {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x7, x3, x2, x5, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x7, x3, x2, x5, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
                  o3data[ix7+x7.size()*(ix3+x3.size()*(ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))] +=
                    -1.0 * fdata[ix6+x4.size()*(ix3)] * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
                }
              }
            }
          }
        }
      }
    }
    out(3)->add_block(o3data, x7, x3, x2, x5, x1, x0);
  }
  {
    std::unique_ptr<double[]> o4data(new double[out(4)->get_size(x7, x6, x2, x5, x1, x0)]);
    std::fill_n(o4data.get(), out(4)->get_size(x7, x6, x2, x5, x1, x0), 0.0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                o4data[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))] +=
                  -1.0 * i0data_sorted[ix7+x7.size()*(ix6+x6.size()*(ix2+x2.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
    }
    out(4)->add_block(o4data, x7, x6, x2, x5, x1, x0);
  }
}

void Task541::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x5 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I704
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x7, x6, c1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
    sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x7.size(), x6.size(), x5.size());
  out()->add_block(odata, x7, x6, x5, x2, x1, x0);
}

void Task542::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I707
  std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, x0, x5, x4, x3)]);
  sort_indices<3,4,0,5,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), x0.size(), x5.size(), x4.size(), x3.size());
  // tensor label (calculated on-the-fly): Gamma254
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix0)] +=
                -1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix4+x2.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x0);
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x2, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x2, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix0)))] +=
                  -1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x2, x0);
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> o1data(new double[out(1)->get_size(x5, x0)]);
      std::fill_n(o1data.get(), out(1)->get_size(x5, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              o1data[ix5+x5.size()*(ix0)] +=
                2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix3+x2.size()*(ix3+x3.size()*(ix4+x1.size()*(ix0)))))];
            }
          }
        }
      }
      out(1)->add_block(o1data, x5, x0);
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x0, x2, x3)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x0, x2, x3), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix0+x0.size()*(ix2+x2.size()*(ix3)))] +=
                  -1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix3+x3.size()*(ix4+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x0, x2, x3);
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x4, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x4, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0)))] +=
                  2.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix3+x2.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x4, x1, x0);
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> o2data(new double[out(2)->get_size(x5, x3, x1, x0)]);
      std::fill_n(o2data.get(), out(2)->get_size(x5, x3, x1, x0), 0.0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                o2data[ix5+x5.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))] +=
                  -1.0 * i0data_sorted[ix5+x5.size()*(ix4+x4.size()*(ix4+x2.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0)))))];
              }
            }
          }
        }
      }
      out(2)->add_block(o2data, x5, x3, x1, x0);
    }
  }
  {
    std::unique_ptr<double[]> o3data(new double[out(3)->get_size(x5, x4, x2, x3, x1, x0)]);
    std::fill_n(o3data.get(), out(3)->get_size(x5, x4, x2, x3, x1, x0), 0.0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data_sorted, o3data, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
    out(3)->add_block(o3data, x5, x4, x2, x3, x1, x0);
  }
}

void Task543::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I707
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    // tensor label: I708
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->add_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task544::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  // tensor label: I708
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1, x0, c2)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1, x0, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,-2,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    dgemm_("T", "N", c2.size(), x2.size()*x1.size()*x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x2, x1, x0, c2);
}

void Task545::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I707
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: I724
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<2,1,0,3,4,5,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task546::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I724
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    dscal_(x5.size()*x4.size()*c1.size()*x3.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,-2,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task547::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I724
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x5, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    dgemm_("T", "N", x3.size(), c1.size()*x5.size()*x4.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x5.size(), x4.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task548::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I707
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
    // tensor label: I848
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<5,3,4,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task549::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I848
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, x0, x1)]);
    sort_indices<1,0,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*x0.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, x1, x0, c1, x2);
}

#endif
