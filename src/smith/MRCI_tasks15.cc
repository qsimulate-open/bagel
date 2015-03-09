//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks15.cc
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

#include <src/smith/MRCI_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task700::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  // tensor label: Gamma32
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I133
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task701::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I134
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task702::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I134
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task703::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I137
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task704::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I137
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,8,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task705::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c5 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c5, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c5, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), a2.size());
        // tensor label: I982
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, c1, a2);
}

void Task706::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I982
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, a4, c3);
    sort_indices<0,1,2,3,1,1,-8,1>(i0data, odata, c1.size(), c5.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c1, c5, a4, c3);
}

void Task707::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c5, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c5, a4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c5.size(), a4.size());
        // tensor label: I985
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c5, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c5, a4, c3)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), c3.size());
        dgemm_("T", "N", a2.size(), c1.size(), c5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c5.size()*a4.size()*c3.size(), i1data_sorted, c5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, c1, a2);
}

void Task708::Task_local::compute() {
  const Index c1 = b(0);
  const Index c5 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I985
  std::unique_ptr<double[]> odata = out()->move_block(c1, c5, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c5, a4, c3);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), c5.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c1, c5, a4, c3);
}

void Task709::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
        // tensor label: I988
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        dgemm_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task710::Task_local::compute() {
  const Index a5 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I988
  std::unique_ptr<double[]> odata = out()->move_block(a5, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,8,1>(i0data, odata, a5.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a5, a2, a4, c3);
}

void Task711::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a5 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a5)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
        // tensor label: I991
        std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a2, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a2, a4, c3)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a2.size(), a4.size(), c3.size());
        dgemm_("T", "N", c1.size(), a2.size(), a5.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, a5.size()*a4.size()*c3.size(), i1data_sorted, a5.size()*a4.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task712::Task_local::compute() {
  const Index a5 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I991
  std::unique_ptr<double[]> odata = out()->move_block(a5, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a5.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a5, a2, a4, c3);
}

void Task713::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, x2, c3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, x2, c3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
        // tensor label: I840
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x3, x2, x1, x0)]);
        sort_indices<2,3,1,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), c3.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, c3.size()*x3.size()*x2.size(), i1data_sorted, c3.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task714::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I840
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma107
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x5, x2, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x5, x2, x4, x1, x0)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I841
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, c3, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, c3, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
      dgemm_("T", "N", x3.size()*x2.size()*x1.size()*x0.size(), c1.size()*c3.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x3.size()*x2.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x3, x2, x1, x0);
}

void Task715::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index c3 = b(2);
  const Index x4 = b(3);
  // tensor label: I841
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, c3, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x5.size(), c3.size(), x4.size());
  }
  out()->put_block(odata, c1, x5, c3, x4);
}

void Task716::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a2, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        // tensor label: I843
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task717::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma278
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x4, x5, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x4, x5, x3, x2, x1, x0)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I844
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        dgemm_("T", "N", x4.size()*x3.size()*x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x4.size()*x3.size()*x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size(), x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x4, x3, x2, x1, x0);
}

void Task718::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I844
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, c1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x7.size(), x6.size(), c1.size(), x5.size());
  }
  out()->put_block(odata, x7, x6, c1, x5);
}

void Task719::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, x2, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x3, x2, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x4.size(), x3.size(), x2.size(), a2.size());
        // tensor label: I846
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, x4, x3, x1, x0)]);
        sort_indices<2,3,1,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x2.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x2.size()*x4.size()*x3.size(), i1data_sorted, x2.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task720::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I846
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x4, x3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x4, x3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x4, x3, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma100
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I847
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x4.size()*x3.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x4.size()*x3.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x2.size(), x4.size(), x3.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x4, x3, x1, x0);
}

void Task721::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I847
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, c1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x7.size(), x6.size(), c1.size(), x5.size());
  }
  out()->put_block(odata, x7, x6, c1, x5);
}

void Task722::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size(), c1.size(), a2.size());
      // tensor label: I849
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task723::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I849
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I850
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x2, x1, x0);
}

void Task724::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c3 = b(2);
  const Index x3 = b(3);
  // tensor label: I850
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c3, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c3, x3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x4.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c3, x3);
}

void Task725::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, c1, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2, c1, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size(), c1.size(), c3.size());
      // tensor label: I852
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task726::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I852
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I853
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x2, x1, x0);
}

void Task727::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c3 = b(2);
  const Index x3 = b(3);
  // tensor label: I853
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c3, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c3, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c3, x3);
}

void Task728::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x5)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
      // tensor label: I855
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task729::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I855
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma221
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x3, x2, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I856
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x1, x0);
}

void Task730::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I856
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), c3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x4, c3, x3, x2);
}

void Task731::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I855
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma104
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I862
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x1, x0);
}

void Task732::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I862
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x4, x3, x2, c3);
}

void Task733::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x5.size());
      // tensor label: I858
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task734::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I858
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma221
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x3, x2, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I859
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x1, x0);
}

void Task735::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I859
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x4.size(), c3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x4, c3, x3, x2);
}

void Task736::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I858
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma104
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I865
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x1, x0);
}

void Task737::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I865
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x4, x3, x2, c3);
}

void Task738::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), x4.size());
        // tensor label: I885
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x5, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x5, x4, x1, x0)]);
        sort_indices<2,1,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task739::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I885
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma240
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I886
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x5, x4, x1, x0);
}

void Task740::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I886
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, c1, c3);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x3.size(), x2.size(), c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3, x3, x2);
}

void Task741::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I885
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x2, x4, x1, x0)]);
      sort_indices<1,2,0,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I892
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, c3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x5, x4, x1, x0);
}

void Task742::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I892
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, c1, x3, x2, c3);
}

void Task743::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I885
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x5, x4, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma296
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x3, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x3, x4, x1, x0)]);
      sort_indices<1,2,0,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x3.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I898
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), c1.size(), x2.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c3.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<5,4,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, c1, c3, x5, x4, x1, x0);
}

void Task744::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I898
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, c1, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c3, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), c3.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x3, c3, c1, x2);
}

void Task745::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma240
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
          sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
          // tensor label: I888
          std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3, x2, x5, c1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3, x2, x5, c1, x4)]);
          sort_indices<3,5,1,2,0,4,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size(), x2.size(), x5.size(), c1.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task746::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c1 = b(4);
  const Index x4 = b(5);
  // tensor label: I888
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3, x2, x5, c1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3, x2, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3, x2, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c1.size(), x4.size());
    // tensor label: I889
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size(), x3.size(), x2.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), a2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, a2, x3, x2, x5, c1, x4);
}

void Task747::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I889
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a3, a2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x3.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2, x3, x2);
}

void Task748::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c1 = b(4);
  const Index x4 = b(5);
  // tensor label: I888
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3, x2, x5, c1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3, x2, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3, x2, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x5.size(), x4.size());
    // tensor label: I919
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size(), x3.size(), x2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), a2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, a2, x3, x2, x5, c1, x4);
}

void Task749::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I919
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a2, a3, x2);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x3.size(), a2.size(), a3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x3, x2, a3, a2);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x3.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2, x3, x2);
}

#endif
