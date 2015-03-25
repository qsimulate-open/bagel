//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks9.cc
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

#include <src/smith/MRCI_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task400::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I537
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma176
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x0, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x0, x4, x2, x1)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x0.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I538
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c1.size());
  out()->put_block(odata, c1, x5, x0, x4);
}

void Task401::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I538
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c1.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, c1, x3, x2, x1);
}

void Task402::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I537
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma178
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x0, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x0, x4, x3, x2)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x0.size(), x4.size(), x3.size(), x2.size());
        // tensor label: I544
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c1.size());
  out()->put_block(odata, c1, x5, x0, x4);
}

void Task403::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I544
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task404::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I540
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task405::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I540
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I541
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x4, x0);
}

void Task406::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I541
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c3.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, c3, x3, x2, x1);
}

void Task407::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I540
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I547
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x4, x0);
}

void Task408::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I547
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c3, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c3, x1);
}

void Task409::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), c4.size());
      // tensor label: I549
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*c1.size(), a2.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task410::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I549
  std::unique_ptr<double[]> odata = out()->move_block(a2, c4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I550
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c4.size());
  out()->put_block(odata, a2, c4, x0, x1);
}

void Task411::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I550
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c4, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c4.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c4, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c4.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c4, x2);
}

void Task412::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c1.size(), a2.size());
      // tensor label: I552
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x1, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task413::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I552
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I553
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c3.size());
  out()->put_block(odata, a4, c3, x1, x0);
}

void Task414::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I553
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), a4.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c3, x2);
}

void Task415::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I552
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I583
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a4.size());
  out()->put_block(odata, a4, c3, x1, x0);
}

void Task416::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I583
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x3, x2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c3.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c3, a4, x3, x2);
}

void Task417::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), c4.size());
      // tensor label: I555
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task418::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I555
  std::unique_ptr<double[]> odata = out()->move_block(a2, c4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I556
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c4.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c4.size());
  out()->put_block(odata, a2, c4, x1, x0);
}

void Task419::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I556
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c4, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c4, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), c4.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c4, x2);
}

void Task420::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I555
  std::unique_ptr<double[]> odata = out()->move_block(a2, c4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I586
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c4.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c4.size(), a2.size());
  out()->put_block(odata, a2, c4, x1, x0);
}

void Task421::Task_local::compute() {
  const Index c4 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I586
  std::unique_ptr<double[]> odata = out()->move_block(c4, a2, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c4.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c4, a2, x3, x2);
}

void Task422::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, a4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), a4.size(), a2.size());
      // tensor label: I558
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task423::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I559
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c3.size());
  out()->put_block(odata, a4, c3, x1, x0);
}

void Task424::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I559
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a4.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c3, x2);
}

void Task425::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I589
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a4.size());
  out()->put_block(odata, a4, c3, x1, x0);
}

void Task426::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I589
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c3, a4, x3, x2);
}

void Task427::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c3.size(), a2.size());
      // tensor label: I561
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1, x1, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task428::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I561
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I562
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c1.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c1.size());
  out()->put_block(odata, a4, c1, x1, x0);
}

void Task429::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I562
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c1, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a4.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c1, x2);
}

void Task430::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I561
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I592
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), a4.size());
  out()->put_block(odata, a4, c1, x1, x0);
}

void Task431::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I592
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c1, a4, x3, x2);
}

void Task432::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, a4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), a4.size(), a2.size());
      // tensor label: I564
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1, x0, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task433::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I564
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I565
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a4.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a4.size(), c1.size());
  out()->put_block(odata, a4, c1, x0, x1);
}

void Task434::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, c1, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a4.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a4, c1, x2);
}

void Task435::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I567
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task436::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I567
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I568
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x5, x4, x0);
}

void Task437::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I568
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, c1, x3, x2, x1);
}

void Task438::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I567
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I574
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x5, x4, x0);
}

void Task439::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I574
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task440::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I570
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task441::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I570
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I571
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x4, x0);
}

void Task442::Task_local::compute() {
  const Index c3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I571
  std::unique_ptr<double[]> odata = out()->move_block(c3, x3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, c3, x3, x2, x1);
}

void Task443::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I570
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I577
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->put_block(odata, c3, x5, x4, x0);
}

void Task444::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I577
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c3, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c3, x1);
}

void Task445::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2, c3, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size(), c3.size(), x1.size());
      // tensor label: I597
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task446::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I597
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma196
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x0, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x0.size(), x1.size());
        // tensor label: I598
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x0.size()*x1.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x0.size(), x1.size(), a2.size());
  out()->put_block(odata, a2, x2, x0, x1);
}

void Task447::Task_local::compute() {
  const Index x5 = b(0);
  const Index a2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I598
  std::unique_ptr<double[]> odata = out()->move_block(x5, a2, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a2.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a2, x4, x3);
}

void Task448::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I642
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task449::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I643
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c1.size(), a4.size());
  out()->put_block(odata, c1, a4, x3, x0);
}

#endif
