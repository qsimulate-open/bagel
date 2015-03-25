//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks8.cc
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

#include <src/smith/MRCI_tasks8.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task350::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I423
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I424
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), c3.size());
  out()->put_block(odata, c3, x0, x2, x1);
}

void Task351::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c3 = b(2);
  const Index x3 = b(3);
  // tensor label: I424
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c3, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c3, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x4.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c3, x3);
}

void Task352::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), c3.size(), a2.size());
      // tensor label: I426
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task353::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I426
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I427
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), c1.size());
  out()->put_block(odata, c1, x0, x2, x1);
}

void Task354::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I427
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c1, x3);
}

void Task355::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size(), c3.size(), a2.size());
    // tensor label: I429
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task356::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I429
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I430
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c4, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c4.size(), x1.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task357::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c4 = b(2);
  const Index x1 = b(3);
  // tensor label: I430
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c4, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c4, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), c4.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c4, x1);
}

void Task358::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, c4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), c4.size());
    // tensor label: I432
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task359::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I432
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I433
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c4, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c4.size(), x1.size());
        dgemm_("T", "N", x0.size(), c4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task360::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c4 = b(2);
  const Index x1 = b(3);
  // tensor label: I433
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c4, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c4, x1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x2.size(), c4.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c4, x1);
}

void Task361::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x3.size());
      // tensor label: I435
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, x0, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), x0.size(), x3.size());
      dgemm_("T", "N", a2.size()*c3.size(), c1.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task362::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I435
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I436
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), c1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), c1.size(), c4.size());
  out()->put_block(odata, c1, c4, x0, x3);
}

void Task363::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I436
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, c4, c1, x1);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x2.size(), c4.size(), c1.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, c1, c4);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x2.size(), x1.size(), c1.size(), c4.size());
  }
  out()->put_block(odata, c1, c4, x2, x1);
}

void Task364::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I435
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I454
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x0.size()*x3.size(), c1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), c1.size(), c4.size());
  out()->put_block(odata, c1, c4, x0, x3);
}

void Task365::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I454
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c1, x2, x1, c4);
}

void Task366::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x3.size());
      // tensor label: I438
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, x0, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), x0.size(), x3.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task367::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I438
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I439
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), c1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), c1.size(), c4.size());
  out()->put_block(odata, c1, c4, x0, x3);
}

void Task368::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I439
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, c1, c4);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x2.size(), x1.size(), c1.size(), c4.size());
  }
  out()->put_block(odata, c1, c4, x2, x1);
}

void Task369::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I438
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I457
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c1.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c1.size(), c4.size());
  out()->put_block(odata, c1, c4, x0, x3);
}

void Task370::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I457
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c1, x2, x1, c4);
}

void Task371::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I438
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma155
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x0, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
      // tensor label: I475
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c4.size()*c1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c4.size(), c1.size());
  out()->put_block(odata, c1, c4, x0, x3);
}

void Task372::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I475
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c1, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x2.size(), c4.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x2, c4, c1, x1);
}

void Task373::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x3.size());
      // tensor label: I441
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x0, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x0.size(), x3.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task374::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I441
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I442
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), c3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), c3.size(), c4.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task375::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I442
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, c3, c4);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x2.size(), x1.size(), c3.size(), c4.size());
  }
  out()->put_block(odata, c3, c4, x2, x1);
}

void Task376::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I441
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I460
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), c4.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task377::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I460
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c3, x2, x1, c4);
}

void Task378::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I441
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma155
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x0, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
      // tensor label: I478
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c3, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c3, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c3.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c4.size()*c3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c4.size(), c3.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task379::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I478
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c3, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c3, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x2.size(), c4.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x2, c4, c3, x1);
}

void Task380::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x3.size());
      // tensor label: I444
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2, x0, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size(), x0.size(), x3.size());
      dgemm_("T", "N", c3.size()*c1.size(), a2.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task381::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I444
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I445
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), a4.size()*a2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), a4.size(), a2.size());
  out()->put_block(odata, a4, a2, x0, x3);
}

void Task382::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I445
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a2.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, a4, a2, x2, x1);
}

void Task383::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I444
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I463
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a2.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a2.size());
  out()->put_block(odata, a4, a2, x0, x3);
}

void Task384::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I463
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, x1, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, x1, a2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), x2.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a4, x2, x1, a2);
}

void Task385::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I444
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma155
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x0, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
      // tensor label: I481
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a2.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a2.size(), a4.size());
  out()->put_block(odata, a4, a2, x0, x3);
}

void Task386::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I481
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, a4, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), a2.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, x2, a2, a4, x1);
}

void Task387::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x3.size());
      // tensor label: I447
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x0, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x0.size(), x3.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task388::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I447
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I448
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), c3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), c3.size(), c4.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task389::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I448
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, c3, c4);
    sort_indices<2,3,0,1,1,1,-1,1>(i1data, odata, x2.size(), x1.size(), c3.size(), c4.size());
  }
  out()->put_block(odata, c3, c4, x2, x1);
}

void Task390::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I447
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I466
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), c4.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task391::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I466
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c3, x2, x1, c4);
}

void Task392::Task_local::compute() {
  const Index c3 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I447
  std::unique_ptr<double[]> odata = out()->move_block(c3, c4, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c4, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c4, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma155
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x0, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x0.size(), x1.size());
      // tensor label: I484
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c3, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c3, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c3.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c4.size()*c3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c4.size(), c3.size());
  out()->put_block(odata, c3, c4, x0, x3);
}

void Task393::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I484
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c3, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c3, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), c4.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x2, c4, c3, x1);
}

void Task394::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: I450
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2, x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2, x0, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size(), x0.size(), x3.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task395::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I450
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2, x0, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma197
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I451
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a2, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a2, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size(), x2.size(), x1.size());
      dgemm_("T", "N", x0.size()*x3.size(), a4.size()*a2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), a4.size(), a2.size());
  out()->put_block(odata, a4, a2, x0, x3);
}

void Task396::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I451
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), a2.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, a2, a4, x1);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x2.size(), a2.size(), a4.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, a4, a2);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x2.size(), x1.size(), a4.size(), a2.size());
  }
  out()->put_block(odata, a4, a2, x2, x1);
}

void Task397::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: I450
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2, x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2, x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2, x0, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I469
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a2.size());
      dgemm_("T", "N", x0.size()*x3.size(), a4.size()*a2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x0.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size(), a4.size(), a2.size());
  out()->put_block(odata, a4, a2, x0, x3);
}

void Task398::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I469
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, x1, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, x1, a2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a4.size(), x2.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a4, x2, x1, a2);
}

void Task399::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), x4.size());
      // tensor label: I537
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", a2.size()*c3.size(), c1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

#endif
