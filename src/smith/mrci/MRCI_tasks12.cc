//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks12.cc
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

#include <src/smith/mrci/MRCI_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task550::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma573
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: I1725
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a2, c1, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a2, c1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->add_block(odata, c1, x1, x0, a2);
}

void Task551::Task_local::compute() {
  const Index x7 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x6 = b(3);
  // tensor label: I1725
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, a2, c1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x7, a2, c1, x6), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x7.size(), a2.size(), c1.size(), x6.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x7, x6);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), x7.size(), x6.size());
  }
  out()->add_block(odata, x7, a2, c1, x6);
}

void Task552::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2, x0, a1), 0.0);
  {
    // tensor label: I144
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2, x0, a1);
}

void Task553::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma48
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I145
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, a1, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, a1, x4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), a1.size(), x4.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task554::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x5, a1, x4);
}

void Task555::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x3, x5, a1, x4);
}

void Task556::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, x5, a1, x4);
}

void Task557::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I148
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x5, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x5, x4)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x5.size(), x4.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task558::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", a1.size()*x5.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, a1.size()*x5.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, a1, x5, x4);
}

void Task559::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, c2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x3, a1, x5, x4);
}

void Task560::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x5.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, a1)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x3, a1, x5, x4);
}

void Task561::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, c2)]);
      sort_indices<3,2,0,1,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x3, a1, x5, x4);
}

void Task562::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<0,3,1,2,0,1,-1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, a1, x5, x4);
}

void Task563::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I151
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task564::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task565::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,2,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task566::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I1075
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task567::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1075
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->add_block(odata, a3, c2, x4, x3);
}

void Task568::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I1078
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task569::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1078
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->add_block(odata, a3, c2, x4, x3);
}

void Task570::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, a2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<1,2,0,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task571::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, a2, a1)]);
      sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, a3.size(), x3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*a2.size(),
             1.0, i0data_sorted, a3.size()*a2.size(), i1data_sorted, a3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task572::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    // tensor label: I154
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size());
    dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task573::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task574::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: h1
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task575::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, c2)]);
        sort_indices<3,2,1,0,0,1,-4,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", a1.size(), x3.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->add_block(odata, x3, a1);
}

void Task576::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, c2)]);
        sort_indices<3,1,2,0,0,1,2,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", a1.size(), x3.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->add_block(odata, x3, a1);
}

void Task577::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<0,3,2,1,0,1,2,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task578::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, x3, a1);
}

void Task579::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a1, x3, c2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, a1, x3, c2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a1.size(), x3.size(), c2.size());
        // tensor label: I1026
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x4, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x4, x0, x2, x1)]);
        sort_indices<2,1,0,3,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x3.size()*x4.size(),
               1.0, i0data_sorted, c2.size()*x3.size()*x4.size(), i1data_sorted, c2.size()*x3.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task580::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1026
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x4, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x4, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x4, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x4, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma339
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x3, x5, x4, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x3, x5, x4, x0, x2, x1)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x3.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c2, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c2, x5)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
        dgemm_("T", "N", x3.size()*x4.size()*x0.size()*x2.size()*x1.size(), c2.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x3.size()*x4.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x3.size(), x4.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x3, x4, x0, x2, x1);
}

void Task581::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma340
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1029
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a1, x5)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x5.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task582::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1029
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, a1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, a1, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, a1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x5)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, c2)]);
      sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), c2.size());
      dgemm_("T", "N", a1.size()*x5.size(), x4.size()*x3.size(), c3.size()*c2.size(),
             1.0, i0data_sorted, c3.size()*c2.size(), i1data_sorted, c3.size()*c2.size(),
             1.0, odata_sorted, a1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, a1, x5);
}

void Task583::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x6 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, c2, x6);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, c2, x6)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), c2.size(), x6.size());
        // tensor label: I1032
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x0, x6, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x0, x6, x2, x1)]);
        sort_indices<1,0,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x0.size(), x6.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x7.size()*x6.size(),
               1.0, i0data_sorted, c2.size()*x7.size()*x6.size(), i1data_sorted, c2.size()*x7.size()*x6.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task584::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1032
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma341
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x5, x6, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x5, x6, x4, x3, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x5.size(), x6.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x0, x6, x2, x1);
}

void Task585::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1032
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma342
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x3, x6, x5, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x3, x6, x5, x4, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x3.size(), x6.size(), x5.size(), x4.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c2)]);
        sort_indices<2,0,1,3,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x0, x6, x2, x1);
}

void Task586::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x7 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x7, x6);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x7, x6)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
        // tensor label: I1044
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x6, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x6, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x6.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x7.size()*x6.size(),
               1.0, i0data_sorted, c2.size()*x7.size()*x6.size(), i1data_sorted, c2.size()*x7.size()*x6.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task587::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1044
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma345
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x0, x4, x3, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task588::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1044
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma346
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,3,4,0,1,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c2)]);
        sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task589::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma349
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
            // tensor label: I1056
            std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,5,1,2,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task590::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1056
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I1057
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, a1, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, a1, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), a1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task591::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1057
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, a1, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a2, a1, x4, x3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a2.size(), a1.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a2, a1);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x4.size(), x3.size(), a2.size(), a1.size());
  }
  out()->add_block(odata, a2, a1, x4, x3);
}

void Task592::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1056
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<1,4,5,0,2,3,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task593::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma350
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4, x6, x5, x3, x0, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x4, x6, x5, x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I1059
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a1, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a1, x7, x6, x5)]);
            sort_indices<3,0,4,5,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task594::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1059
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, a1, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, a1, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, a1, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x4, x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x4, x3, a1)]);
    sort_indices<0,1,2,3,0,1,1,2>(i1data, i1data_sorted, a2.size(), x4.size(), x3.size(), a1.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*x3.size()*a1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, x4, x3, a1, x7, x6, x5);
}

void Task595::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma351
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x5, x4, x0, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x6, x5, x4, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I1062
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task596::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1062
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, a1, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, a1, x3, x7, x6, x5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, a1, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, a2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,2>(i1data, i1data_sorted, x4.size(), a1.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*a1.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), a1.size(), x3.size());
  out()->add_block(odata, x4, a1, x3, x7, x6, x5);
}

void Task597::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma359
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x4, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x4, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1086
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x5, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x5, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x5.size(), a1.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task598::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x5 = b(2);
  const Index a1 = b(3);
  // tensor label: I1086
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), a1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, x5, a1);
}

void Task599::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma366
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x3, x5, x4, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x3, x5, x4, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x3.size(), x5.size(), x4.size(), x2.size(), x1.size());
            // tensor label: I1107
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, a1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, a1, x6)]);
            sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), a1.size(), x6.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

#endif
