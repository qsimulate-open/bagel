//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks22.cc
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

#include <src/smith/MRCI_tasks22.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1050::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I213
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1051::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I213
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I214
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1052::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I214
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task1053::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I213
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1465
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1054::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1465
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), c4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, c4, x2, x1);
}

void Task1055::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I213
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma29
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1471
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, c4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), c4.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1056::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1471
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x3, x2, x1, c4);
}

void Task1057::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I216
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1058::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I217
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1059::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I217
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task1060::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1468
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1061::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1468
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), c4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, c4, x2, x1);
}

void Task1062::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma29
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I1474
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, c4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), c4.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task1063::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1474
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, x3, x2, x1, c4);
}

void Task1064::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I219
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, a3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, a3, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), a3.size(), a1.size());
    dgemm_("T", "N", x0.size(), c2.size()*a3.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1065::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I220
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a3.size()*a1.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1066::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I220
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task1067::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I223
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a1.size()*a3.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*a3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), a3.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1068::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I223
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task1069::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I226
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*c2.size()*a3.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a3.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1070::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I226
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task1071::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I229
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*a3.size()*c2.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*c2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1072::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I229
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task1073::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I232
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*c2.size()*a1.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1074::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  // tensor label: I232
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3);
}

void Task1075::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I235
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*a1.size()*c2.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*c2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1076::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  // tensor label: I235
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3);
}

void Task1077::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c5, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c5, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c5.size(), a3.size());
      // tensor label: I1483
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, a4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, a4, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), a4.size(), a1.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*a1.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1078::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a1 = b(3);
  // tensor label: I1483
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, a4, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, a4, a1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), c5.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, x1, c5, a4, a1);
}

void Task1079::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c5, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c5.size(), a4.size());
      // tensor label: I1486
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, a4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, a4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), a4.size(), a1.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*a1.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1080::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a1 = b(3);
  // tensor label: I1486
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, a4, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, a4, a1);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x1.size(), c5.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, x1, c5, a4, a1);
}

void Task1081::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c5, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c5, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c5.size(), a1.size());
      // tensor label: I1489
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, a4, a3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), a4.size(), a3.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*a3.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1082::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a3 = b(3);
  // tensor label: I1489
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, a4, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, a4, a3);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x1.size(), c5.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, x1, c5, a4, a3);
}

void Task1083::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c5, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c5.size(), a4.size());
      // tensor label: I1492
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, a4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), a4.size(), a3.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*a3.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1084::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a3 = b(3);
  // tensor label: I1492
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, a4, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, a4, a3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), c5.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, x1, c5, a4, a3);
}

void Task1085::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a5, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a5, c4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a5.size(), c4.size(), a3.size());
      // tensor label: I1495
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), a5.size(), c4.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*a1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1086::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1495
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, a5, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x1.size(), a1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, x1, a1, a5, c4);
}

void Task1087::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a5.size());
      // tensor label: I1498
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), a5.size(), c4.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*a1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1088::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1498
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, a5, c4);
    sort_indices<0,1,2,3,1,1,8,1>(i0data, odata, x1.size(), a1.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, x1, a1, a5, c4);
}

void Task1089::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a5, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a5, c4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a5.size(), c4.size(), a1.size());
      // tensor label: I1501
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a3, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a3, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), a5.size(), c4.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*a3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1090::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1501
  std::unique_ptr<double[]> odata = out()->move_block(x1, a3, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, a5, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), a3.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, x1, a3, a5, c4);
}

void Task1091::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a5.size());
      // tensor label: I1504
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a3, a5, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a3, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), a3.size(), a5.size(), c4.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*a3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1092::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  const Index a5 = b(2);
  const Index c4 = b(3);
  // tensor label: I1504
  std::unique_ptr<double[]> odata = out()->move_block(x1, a3, a5, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, a5, c4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x1.size(), a3.size(), a5.size(), c4.size());
  }
  out()->put_block(odata, x1, a3, a5, c4);
}

void Task1093::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c5, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c5, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c5.size(), a3.size());
      // tensor label: I1579
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), a1.size());
      dgemm_("T", "N", x1.size()*a3.size(), c2.size()*a1.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*a3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1094::Task_local::compute() {
  const Index c2 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a1 = b(3);
  // tensor label: I1579
  std::unique_ptr<double[]> odata = out()->move_block(c2, c5, a4, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c5, a4, a1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), c5.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, c2, c5, a4, a1);
}

void Task1095::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c5, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c5.size(), a4.size());
      // tensor label: I1582
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), a1.size());
      dgemm_("T", "N", x1.size()*a3.size(), c2.size()*a1.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*a3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1096::Task_local::compute() {
  const Index c2 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a1 = b(3);
  // tensor label: I1582
  std::unique_ptr<double[]> odata = out()->move_block(c2, c5, a4, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c5, a4, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), c5.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, c2, c5, a4, a1);
}

void Task1097::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c5, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c5, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c5.size(), a1.size());
      // tensor label: I1585
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, a3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), a3.size());
      dgemm_("T", "N", x1.size()*a1.size(), c2.size()*a3.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task1098::Task_local::compute() {
  const Index c2 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index a3 = b(3);
  // tensor label: I1585
  std::unique_ptr<double[]> odata = out()->move_block(c2, c5, a4, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c5, a4, a3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), c5.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, c2, c5, a4, a3);
}

void Task1099::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I219
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c5, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c5.size(), a4.size());
      // tensor label: I1588
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c5, a4, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c5, a4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c5.size(), a4.size(), a3.size());
      dgemm_("T", "N", x1.size()*a1.size(), c2.size()*a3.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, x1.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

#endif
