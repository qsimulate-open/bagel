//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks13.cc
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

#include <src/smith/MRCI_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task600::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma244
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x1, x2)]);
          sort_indices<0,1,2,5,3,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
          // tensor label: I744
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, c2, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, c2, x5, x4)]);
          sort_indices<4,5,0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), c2.size(), x5.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task601::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index c2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I744
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, c2, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x2, c2, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x2, c2, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
    // tensor label: I745
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), x2.size());
    dgemm_("T", "N", c2.size()*x5.size()*x4.size(), x3.size()*a1.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c2.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x5.size(), x4.size(), x3.size(), a1.size(), x2.size());
  out()->put_block(odata, x3, a1, x2, c2, x5, x4);
}

void Task602::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I745
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, a3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), a1.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, a1, a3, x2);
}

void Task603::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I759
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x0, x6, x5, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x0, x6, x5, x1)]);
        sort_indices<1,3,4,0,2,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x0.size(), x6.size(), x5.size(), x1.size());
        dgemm_("T", "N", a1.size(), c2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task604::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I759
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x0, x6, x5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma250
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x1, x4, x3, x2)]);
        sort_indices<5,6,7,0,1,2,3,4,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: I760
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x0, x6, x5, x1);
}

void Task605::Task_local::compute() {
  const Index c2 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I760
  std::unique_ptr<double[]> odata = out()->move_block(c2, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c2, x4, x3, x2);
}

void Task606::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I759
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x0, x6, x5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma251
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x1, x2)]);
        sort_indices<4,5,7,0,1,2,3,6,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I763
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c2, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c2, x2)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c2.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x0, x6, x5, x1);
}

void Task607::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I763
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c2, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, c2, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), x3.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, c2, x2);
}

void Task608::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I765
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task609::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I765
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I766
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x2, x1, x0);
}

void Task610::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I766
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task611::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, a3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, a3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), a3.size(), a1.size());
      // tensor label: I768
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x0.size()*x1.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), x1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task612::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I768
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I769
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a3.size());
  out()->put_block(odata, a3, x0, x1, x2);
}

void Task613::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I769
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task614::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I807
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task615::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I807
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I808
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task616::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I808
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, a3, x4, x3, x2);
}

void Task617::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I807
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I814
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task618::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I814
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, a3, x2);
}

void Task619::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I810
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x0.size(), x1.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task620::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I810
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma230
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: I811
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), a3.size());
  out()->put_block(odata, a3, x5, x0, x1);
}

void Task621::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I811
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, a3, x4, x3, x2);
}

void Task622::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I810
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I817
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), a3.size());
  out()->put_block(odata, a3, x5, x0, x1);
}

void Task623::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I817
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, a3, x2);
}

void Task624::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma276
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x1, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x2, x1, x3)]);
          sort_indices<0,2,3,5,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x2.size(), x1.size(), x3.size());
          // tensor label: I837
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, x5, a1, x4)]);
          sort_indices<3,5,2,1,0,4,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task625::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I837
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, x2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I838
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), c2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->put_block(odata, c2, x3, x2, x5, a1, x4);
}

void Task626::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I838
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c2.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, c2, x3, a3, x2);
}

void Task627::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma568
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x1.size(), x4.size());
      // tensor label: I1717
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task628::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x4 = b(3);
  // tensor label: I1717
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, c2, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a1.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, x5, a1, c2, x4);
}

void Task629::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma569
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I1719
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task630::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1719
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x5, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), a1.size(), x5.size(), x4.size());
  }
  out()->put_block(odata, c2, a1, x5, x4);
}

void Task631::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma572
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x1, x6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x1.size(), x6.size());
      // tensor label: I1725
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, c2, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, c2, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), c2.size(), x6.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task632::Task_local::compute() {
  const Index x7 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x6 = b(3);
  // tensor label: I1725
  std::unique_ptr<double[]> odata = out()->move_block(x7, a1, c2, x6);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, c2, x6);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x7.size(), a1.size(), c2.size(), x6.size());
  }
  out()->put_block(odata, x7, a1, c2, x6);
}

void Task633::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma573
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: I1727
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x7, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x7, x6)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task634::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x7 = b(2);
  const Index x6 = b(3);
  // tensor label: I1727
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x7, x6);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x7, x6);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), a1.size(), x7.size(), x6.size());
  }
  out()->put_block(odata, c2, a1, x7, x6);
}

void Task635::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a2);
  {
    // tensor label: I108
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, x0, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), x1.size(), x0.size(), a2.size());
  }
  out()->put_block(odata, x0, x1, c1, a2);
}

void Task636::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I109
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task637::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I109
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I110
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
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task638::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I110
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c1, x3);
}

void Task639::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I112
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, c1, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), c1.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task640::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: I113
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", a2.size()*c1.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task641::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I113
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task642::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I116
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", c1.size()*a2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x3.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task643::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I116
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task644::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x3.size());
      // tensor label: I868
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      dgemm_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task645::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I868
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c1, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c1, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), c4.size(), c1.size(), c3.size());
  }
  out()->put_block(odata, x2, c4, c1, c3);
}

void Task646::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x3.size());
      // tensor label: I871
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c1, c3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c1.size(), c3.size());
      dgemm_("T", "N", a2.size()*x3.size(), x2.size()*c1.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task647::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I871
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c1, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c1, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x2.size(), c4.size(), c1.size(), c3.size());
  }
  out()->put_block(odata, x2, c4, c1, c3);
}

void Task648::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c1.size(), x3.size());
      // tensor label: I874
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task649::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I874
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, x2, c4, a3, a2);
}

#endif
