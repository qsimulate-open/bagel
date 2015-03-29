//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks27.cc
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

#include <src/smith/MRCI_tasks27.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1300::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), a2.size());
      // tensor label: I1644
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x0.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), c3.size()*x5.size(),
             1.0, i0data_sorted, c3.size()*x5.size(), i1data_sorted, c3.size()*x5.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1301::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1644
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma526
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: I1645
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), c3.size());
  out()->put_block(odata, c3, x5, x0, x1);
}

void Task1302::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1645
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), c3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x4, c3, x3, x2);
}

void Task1303::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1644
  std::unique_ptr<double[]> odata = out()->move_block(c3, x5, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x0, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1648
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, c3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), c3.size());
        dgemm_("T", "N", x5.size()*x0.size()*x1.size(), c3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x1.size(), c3.size());
  out()->put_block(odata, c3, x5, x0, x1);
}

void Task1304::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1648
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x4, x3, x2, c3);
}

void Task1305::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1650
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, x3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, x3, a1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), x3.size(), a1.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*a1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1306::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index a1 = b(3);
  // tensor label: I1650
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), a1.size());
      // tensor label: I1651
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, a2, x3, a1);
}

void Task1307::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1651
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, x2, c4, a3, a2);
}

void Task1308::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma526
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
          sort_indices<0,2,4,5,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
          // tensor label: I1662
          std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3, x2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3, x2, x5, a1, x4)]);
          sort_indices<3,5,1,2,0,4,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size(), x2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1309::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1662
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3, x2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I1663
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size(), x3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), a2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, a2, x3, x2, x5, a1, x4);
}

void Task1310::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1663
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a3, a2);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x3.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2, x3, x2);
}

void Task1311::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma50
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
          // tensor label: I1665
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a2, x5, a1, x4)]);
          sort_indices<3,5,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1312::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1665
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I1666
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x3, x2, a2, x5, a1, x4);
}

void Task1313::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I1666
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, x2, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), x3.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a3, x3, x2, a2);
}

void Task1314::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma545
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x2, x3, x1)]);
          sort_indices<0,2,3,4,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
          // tensor label: I1668
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, x5, a1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, x5, a1, x4)]);
          sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x5.size(), a1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), a2.size()*a1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1315::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  const Index x4 = b(5);
  // tensor label: I1668
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, x2, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2, x2, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, x2, x5, a1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
    // tensor label: I1669
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size()*a2.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size(), a2.size(), x2.size());
  out()->put_block(odata, x3, a2, x2, x5, a1, x4);
}

void Task1316::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1669
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, a3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, a3, x2);
}

void Task1317::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, c1, x0);
  {
    // tensor label: I260
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x1, x0);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c2, x1, c1, x0);
}

void Task1318::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I260
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I261
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x1, x0);
}

void Task1319::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x3, c4, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x3, c4, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x3.size(), c4.size(), x2.size());
      // tensor label: I262
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, c2, c3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), c2.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*c2.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x3, x2);
}

void Task1320::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I262
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, c2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, c2, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), c4.size(), c2.size(), c3.size());
  }
  out()->put_block(odata, c1, c4, c2, c3);
}

void Task1321::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I261
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), a4.size());
      // tensor label: I298
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, a3, x2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), a3.size(), x2.size());
      dgemm_("T", "N", c1.size()*c2.size(), x3.size()*x2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, c2, x3, x2);
}

void Task1322::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I298
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, a4, x3, a3, x2);
}

void Task1323::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I260
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma548
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5, x1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5, x1, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x1.size(), x4.size());
      // tensor label: I1677
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, c2, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, c2, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x1, x0);
}

void Task1324::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index c2 = b(2);
  const Index x4 = b(3);
  // tensor label: I1677
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, c2, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x5.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, c1, x5, c2, x4);
}

void Task1325::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I260
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma549
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x7, x1, x6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x7, x1, x6)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x7.size(), x1.size(), x6.size());
      // tensor label: I1679
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x7, c2, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x7, c2, x6)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x7.size(), c2.size(), x6.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x1, x0);
}

void Task1326::Task_local::compute() {
  const Index c1 = b(0);
  const Index x7 = b(1);
  const Index c2 = b(2);
  const Index x6 = b(3);
  // tensor label: I1679
  std::unique_ptr<double[]> odata = out()->move_block(c1, x7, c2, x6);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x7, c2, x6);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), x7.size(), c2.size(), x6.size());
  }
  out()->put_block(odata, c1, x7, c2, x6);
}

void Task1327::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I1112
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3, a2, a4);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c3.size(), a2.size(), a4.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task1328::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, x0, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, x0, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), x0.size(), a4.size());
      // tensor label: I1113
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1329::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1113
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I1114
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

void Task1330::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1114
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x3.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, c1, x3, c3, x2);
}

void Task1331::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c6, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c6, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c6.size(), a2.size());
      // tensor label: I1290
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c6, c3, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c6, c3, c5)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c6.size(), c3.size(), c5.size());
      dgemm_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), c6.size()*c5.size(),
             1.0, i0data_sorted, c6.size()*c5.size(), i1data_sorted, c6.size()*c5.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1332::Task_local::compute() {
  const Index c1 = b(0);
  const Index c6 = b(1);
  const Index c3 = b(2);
  const Index c5 = b(3);
  // tensor label: I1290
  std::unique_ptr<double[]> odata = out()->move_block(c1, c6, c3, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c6, c3, c5);
    sort_indices<0,1,2,3,1,1,8,1>(i0data, odata, c1.size(), c6.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c1, c6, c3, c5);
}

void Task1333::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c6 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a2, c6, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a2, c6, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c6.size(), a4.size());
      // tensor label: I1292
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c6, c3, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c6, c3, c5)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c6.size(), c3.size(), c5.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), c6.size()*c5.size(),
             1.0, i0data_sorted, c6.size()*c5.size(), i1data_sorted, c6.size()*c5.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1334::Task_local::compute() {
  const Index c1 = b(0);
  const Index c6 = b(1);
  const Index c3 = b(2);
  const Index c5 = b(3);
  // tensor label: I1292
  std::unique_ptr<double[]> odata = out()->move_block(c1, c6, c3, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c6, c3, c5);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), c6.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c1, c6, c3, c5);
}

void Task1335::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& a6 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a6, c3, a5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a6, c3, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a6.size(), c3.size(), a5.size());
      // tensor label: I1310
      std::unique_ptr<double[]> i1data = in(1)->get_block(a6, a2, a5, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a6, a2, a5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a6.size(), a2.size(), a5.size(), a4.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*a4.size(), a6.size()*a5.size(),
             1.0, i0data_sorted, a6.size()*a5.size(), i1data_sorted, a6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1336::Task_local::compute() {
  const Index a6 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index a4 = b(3);
  // tensor label: I1310
  std::unique_ptr<double[]> odata = out()->move_block(a6, a2, a5, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a6, a2, a5, a4);
    sort_indices<0,1,2,3,1,1,8,1>(i0data, odata, a6.size(), a2.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a6, a2, a5, a4);
}

void Task1337::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a6 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a6);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a6)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a6.size());
      // tensor label: I1312
      std::unique_ptr<double[]> i1data = in(1)->get_block(a6, a2, a5, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a6, a2, a5, a4)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a6.size(), a2.size(), a5.size(), a4.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*a4.size(), a6.size()*a5.size(),
             1.0, i0data_sorted, a6.size()*a5.size(), i1data_sorted, a6.size()*a5.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1338::Task_local::compute() {
  const Index a6 = b(0);
  const Index a2 = b(1);
  const Index a5 = b(2);
  const Index a4 = b(3);
  // tensor label: I1312
  std::unique_ptr<double[]> odata = out()->move_block(a6, a2, a5, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a6, a2, a5, a4);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, a6.size(), a2.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a6, a2, a5, a4);
}

void Task1339::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, x2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), x2.size(), a4.size());
      // tensor label: I1356
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1340::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1356
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1357
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x1, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x1, c3, x0)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x1.size(), c3.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*c3.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x3, x2);
}

void Task1341::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1357
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, c3, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), x1.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, c1, x1, c3, x0);
}

void Task1342::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  // scalar
  // tensor label: Gamma558
  std::unique_ptr<double[]> i0data = in(0)->get_block();
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I1697
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  dgemm_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1343::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1697
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task1344::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1112
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, a2, a4), 0.0);
  // scalar
  // tensor label: Gamma560
  std::unique_ptr<double[]> i0data = in(0)->get_block();
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size()]);
  sort_indices<0,1,1,1>(i0data, i0data_sorted);
  // tensor label: I1701
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  dgemm_("T", "N", 1, c1.size()*a4.size()*c3.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, 1);
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, c3, a2, a4);
}

void Task1345::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1701
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task1346::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I1640
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, a1, a2);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a1.size(), a2.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task1347::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1640
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c4.size(), a2.size());
      // tensor label: I1641
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, c3, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", a1.size()*a2.size(), x1.size()*x0.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a1, a2);
}

void Task1348::Task_local::compute() {
  const Index c4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1641
  std::unique_ptr<double[]> odata = out()->move_block(c4, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1642
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, x2, c3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), x2.size(), c3.size());
      dgemm_("T", "N", x1.size()*x0.size(), c4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c4.size(), c3.size());
  out()->put_block(odata, c4, c3, x1, x0);
}

void Task1349::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1642
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, x2, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), c4.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x3, c4, x2, c3);
}

#endif
