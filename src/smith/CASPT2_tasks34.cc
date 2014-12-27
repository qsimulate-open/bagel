//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks34.cc
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


#include <src/smith/CASPT2_tasks34.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1650::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1775
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, x2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, x2, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), x2.size(), a1.size());
        // tensor label: I2090
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
        sort_indices<2,5,6,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}

void Task1651::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I2090
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma429
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x2, x4, x3, x1, x0);
}

void Task1652::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
          // tensor label: I1783
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x6, x5, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x6, x5, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x6.size(), x5.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x7.size()*x6.size()*x5.size()*a1.size(),
                 1.0, i0data_sorted, x7.size()*x6.size()*x5.size()*a1.size(), i1data_sorted, x7.size()*x6.size()*x5.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1653::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index a1 = b(4);
  // tensor label: I1783
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x6, x5, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x7, x6, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x7, x6, x5, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        // tensor label: I1784
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x7, x0, x6, x5, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x7.size()*x6.size()*x5.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x7.size(), x6.size(), x5.size());
  out()->put_block(odata, ci0, x7, x6, x5, a1);
}

void Task1654::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1784
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  {
    // tensor label: Gamma436
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}

void Task1655::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
          // tensor label: I1790
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x3.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1656::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1790
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I1791
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
      sort_indices<5,4,0,1,2,3,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*x1.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a1);
}

void Task1657::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1791
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1792
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a3, c2, a1);
}

void Task1658::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1792
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1659::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1791
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1796
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,5,6,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a3, c2, a1);
}

void Task1660::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1796
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1661::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1790
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size());
    // tensor label: I2123
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a1);
}

void Task1662::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I2123
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1663::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
          // tensor label: I1840
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, a3, c2, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), a3.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*a3.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x3.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x3.size()*a3.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1664::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I1841
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a3, c2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x3.size()*a3.size()*c2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x3.size(), a3.size(), c2.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1665::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1841
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a3, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), x1.size());
      // tensor label: I1842
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a3, c2);
}

void Task1666::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1842
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1667::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1841
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a3, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x0.size(), x1.size());
      // tensor label: I1850
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a3, c2);
}

void Task1668::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1850
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1669::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a3.size());
    // tensor label: I1845
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, c2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), c2.size());
    dgemm_("T", "N", a3.size(), ci0.size()*x3.size()*a1.size()*c2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x3.size(), a1.size(), c2.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1670::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1845
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1846
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1671::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1846
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1672::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1845
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1854
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1673::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1854
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1674::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I1857
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, a3);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, a3)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), a3.size());
  dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*a3.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), a3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1675::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  // tensor label: I1857
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), x1.size(), x2.size());
        // tensor label: I1858
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
        sort_indices<3,4,2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a3.size(), ci0.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,2,0,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a3);
}

void Task1676::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1858
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1677::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I1861
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, a1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), a1.size());
  dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x3.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<2,3,1,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x3.size(), a1.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1678::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  // tensor label: I1861
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        // tensor label: I1862
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
        sort_indices<3,4,2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a1);
}

void Task1679::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1862
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1680::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1873
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1681::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  // tensor label: I1873
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma459
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task1682::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1876
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1683::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  // tensor label: I1876
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma459
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task1684::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1840
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I1903
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x3.size()*a1.size()*a3.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x3.size(), a1.size(), a3.size());
  out()->put_block(odata, ci0, x3, a3, c2, a1);
}

void Task1685::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);
  // tensor label: I1903
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());
      // tensor label: I1904
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a3);
}

void Task1686::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1904
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1687::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
          // tensor label: I1864
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c2, a3, a1)]);
          sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x1.size()*c2.size()*a3.size()*a1.size(),
                 1.0, i0data_sorted, x1.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x1.size()*c2.size()*a3.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1688::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, c4)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), c4.size());
      // tensor label: I1865
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
      sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*c2.size()*a3.size()*a1.size(), x0.size()*c4.size(),
             1.0, i0data_sorted, x0.size()*c4.size(), i1data_sorted, x0.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1689::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);
  // tensor label: I1865
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
  // tensor label: I1866
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}

void Task1690::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1866
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1691::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);
  // tensor label: I1865
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
  // tensor label: I1870
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}

void Task1692::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1870
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1693::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());
    // tensor label: I1879
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a3, c4, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a3.size(), c4.size(), a1.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x1.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x1.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1694::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);
  // tensor label: I1879
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I1880
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a3, c4, a1);
}

void Task1695::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1880
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1696::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);
  // tensor label: I1879
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I1884
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a3, c4, a1);
}

void Task1697::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1884
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1698::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I1887
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a3.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x1.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x1.size(), c2.size(), a3.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1699::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  // tensor label: I1887
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I1888
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a3);
}

