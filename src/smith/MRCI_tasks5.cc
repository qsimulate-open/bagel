//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks5.cc
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

#include <src/smith/MRCI_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task200::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I309
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x7, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x7, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x7.size(), x6.size());
    // tensor label: I364
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<0,4,5,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task201::Task_local::compute() {
  const Index a2 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(a2, x5, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x5, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a2.size(), x5.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, a2, x5, x4, x3);
}

void Task202::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma101
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x3, x5, x2, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x3, x5, x2, x4, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I312
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,2,5,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task203::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I312
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c2, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
    // tensor label: I313
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, c2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), c2.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->put_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task204::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I313
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, c1, x4, x3, c2);
}

void Task205::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma102
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x4, x5, x2, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x4, x5, x2, x3, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I315
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c1, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c1, x3, x7, x6, x5)]);
            sort_indices<3,4,0,5,2,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), c1.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task206::Task_local::compute() {
  const Index x4 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I315
  std::unique_ptr<double[]> odata = out()->move_block(x4, c1, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, c1, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, c1, x3, x7, x6, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c2, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
    // tensor label: I316
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, c1, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c2.size(), c1.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*c1.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), c1.size(), x3.size());
  out()->put_block(odata, x4, c1, x3, x7, x6, x5);
}

void Task207::Task_local::compute() {
  const Index x4 = b(0);
  const Index c2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I316
  std::unique_ptr<double[]> odata = out()->move_block(x4, c2, c1, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c2, c1, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), c2.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x4, c2, c1, x3);
}

void Task208::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma104
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I321
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task209::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x5.size());
      // tensor label: I322
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, c1, x5);
}

void Task210::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I322
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x4, x3);
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
  out()->put_block(odata, a3, c2, x4, x3);
}

void Task211::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x5.size());
      // tensor label: I325
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, c1, x5);
}

void Task212::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I325
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x4, x3);
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
  out()->put_block(odata, a3, c2, x4, x3);
}

void Task213::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I321
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
      // tensor label: I334
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, c1, x5);
}

void Task214::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I334
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, a2, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), c3.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x4, c3, a2, x3);
}

void Task215::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma107
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x5, x2, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x5, x2, x4, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I330
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task216::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I330
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x5.size());
      // tensor label: I331
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, c1, x5);
}

void Task217::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I331
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, a3, x4, x3, c2);
}

void Task218::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma109
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I336
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x5)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task219::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I336
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x5.size());
      // tensor label: I337
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, c1, x5);
}

void Task220::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I337
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, a2, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), c3.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x4, c3, a2, x3);
}

void Task221::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma114
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5, x2, x6, x4, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x5, x2, x6, x4, x3, x1, x0)]);
            sort_indices<0,1,3,4,5,2,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x5.size(), x2.size(), x6.size(), x4.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I351
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, c1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, c1, x6)]);
            sort_indices<3,0,5,1,2,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), c1.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task222::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index c1 = b(4);
  const Index x6 = b(5);
  // tensor label: I351
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x7, c1, x6);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, c1, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
    // tensor label: I352
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*c1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*c1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), c1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x5, x4, x3, x7, c1, x6);
}

void Task223::Task_local::compute() {
  const Index a2 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I352
  std::unique_ptr<double[]> odata = out()->move_block(a2, x5, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x5, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a2.size(), x5.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, a2, x5, x4, x3);
}

void Task224::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma115
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x2, x6, x5, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x2, x6, x5, x4, x1, x0)]);
            sort_indices<0,1,3,4,5,2,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x2.size(), x6.size(), x5.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I354
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, c1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, c1, x6)]);
            sort_indices<3,2,5,0,1,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), c1.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task225::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index c1 = b(4);
  const Index x6 = b(5);
  // tensor label: I354
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x7, c1, x6);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, c1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, c1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, c1, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
    // tensor label: I355
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*c1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*c1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), c1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x5, x4, x3, x7, c1, x6);
}

void Task226::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I355
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, a2, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), x4.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, a2, x3);
}

void Task227::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma119
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x4, x2, x3, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x4, x2, x3, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I366
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c1, x7, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c1, x7, x6)]);
            sort_indices<4,5,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c1.size(), x7.size(), x6.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task228::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  const Index x7 = b(4);
  const Index x6 = b(5);
  // tensor label: I366
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, c1, x7, x6);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x7, x6);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x7, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x7.size(), x6.size());
    // tensor label: I367
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x5, x4, x3, c1, x7, x6);
}

void Task229::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I367
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, a2, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x4.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, a2, x3);
}

void Task230::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma122
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x5, x2, x4, x1, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x6, x5, x2, x4, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I375
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task231::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I375
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I376
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->put_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task232::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I376
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x4, a2, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), x4.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, c1, x4, a2, x3);
}

void Task233::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma550
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        // tensor label: I1681
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task234::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1681
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, c1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), c1.size(), x5.size());
  }
  out()->put_block(odata, x7, x6, c1, x5);
}

void Task235::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma551
        std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x8, x2, x7, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x9, x8, x2, x7, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x8.size(), x2.size(), x7.size(), x1.size(), x0.size());
        // tensor label: I1683
        std::unique_ptr<double[]> i1data = in(1)->get_block(x9, x8, c1, x7);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x9, x8, c1, x7)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x9.size(), x8.size(), c1.size(), x7.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task236::Task_local::compute() {
  const Index x9 = b(0);
  const Index x8 = b(1);
  const Index c1 = b(2);
  const Index x7 = b(3);
  // tensor label: I1683
  std::unique_ptr<double[]> odata = out()->move_block(x9, x8, c1, x7);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x8, c1, x7);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x9.size(), x8.size(), c1.size(), x7.size());
  }
  out()->put_block(odata, x9, x8, c1, x7);
}

void Task237::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0, c1, a2);
  {
    // tensor label: I27
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3, x0, a2);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c3.size(), x0.size(), a2.size());
  }
  out()->put_block(odata, c3, x0, c1, a2);
}

void Task238::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I28
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size(), c1.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task239::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I28
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I29
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x1, x0);
}

void Task240::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I29
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x3.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, c1, x3, c3, x2);
}

void Task241::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I31
  std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task242::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I32
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->put_block(odata, c3, x0);
}

void Task243::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I32
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c3, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c3, x1);
}

void Task244::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I34
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
  dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task245::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I34
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I35
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->put_block(odata, c1, x0);
}

void Task246::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I35
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task247::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I37
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task248::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I37
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I38
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size());
    dgemm_("T", "N", a2.size()*c3.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c3.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), c1.size());
  out()->put_block(odata, c1, a2, c3, x1);
}

void Task249::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  // tensor label: I38
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c1.size(), c4.size());
  }
  out()->put_block(odata, c1, c4);
}

#endif
