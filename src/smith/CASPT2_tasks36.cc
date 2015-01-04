//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks36.cc
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


#include <src/smith/CASPT2_tasks36.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1750::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1876
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1751::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma470
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x4, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
          // tensor label: I1552
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x5, x4)]);
          sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x5.size(), x4.size());
          dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*x5.size()*x4.size(), i1data_sorted, x1.size()*x0.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1752::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1552
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      // tensor label: I1553
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
      dgemm_("T", "N", x5.size()*x4.size(), x1.size()*x0.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1753::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1553
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task1754::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1552
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x5, x4), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I1915
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x4, a1, x5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x4, a1, x5)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x4.size(), a1.size(), x5.size());
      dgemm_("T", "N", x0.size()*x1.size(), x4.size()*x5.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x4.size(), x5.size());
  out()->put_block(odata, x1, x0, x5, x4);
}

void Task1755::Task_local::compute() {
  const Index a2 = b(0);
  const Index x4 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1915
  std::unique_ptr<double[]> odata = out()->move_block(a2, x4, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), a2.size());
  }
  out()->put_block(odata, a2, x4, a1, x5);
}

void Task1756::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma591
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x5, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x2, x5, x4, x3, x1, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I1996
              std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, x5, x4, x3)]);
              sort_indices<0,3,4,5,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), x5.size(), x4.size(), x3.size());
              dgemm_("T", "N", ci0.size(), 1, x2.size()*x1.size()*x0.size()*x5.size()*x4.size()*x3.size(),
                     1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*x5.size()*x4.size()*x3.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*x5.size()*x4.size()*x3.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1757::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1996
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
    // tensor label: I1997
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1758::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1997
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task1759::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma609
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x0, x3, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x0, x3, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
              // tensor label: I2050
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, x5, x0, x1, x2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, x5, x0, x1, x2)]);
              sort_indices<2,1,3,0,5,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), x5.size(), x0.size(), x1.size(), x2.size());
              dgemm_("T", "N", ci0.size(), 1, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(), i1data_sorted, x3.size()*x4.size()*x5.size()*x0.size()*x1.size()*x2.size(),
                     1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1760::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I2050
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4, x5, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x4, x5, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x4, x5, x0, x1, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, x1, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), x1.size(), x2.size());
    // tensor label: I2051
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c1, x4, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c1, x4, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c1.size(), x4.size(), x5.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
  out()->put_block(odata, x3, x4, x5, x0, x1, x2);
}

void Task1761::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index x4 = b(2);
  const Index x5 = b(3);
  // tensor label: I2051
  std::unique_ptr<double[]> odata = out()->move_block(x3, c1, x4, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x3, c1, x4, x5);
}

