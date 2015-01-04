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
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1419
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task1651::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1417
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, x4, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, x4, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, x4, x2, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), x1.size(), x2.size());
    // tensor label: I1756
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x5, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x5, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x5.size(), a2.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<5,4,3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->put_block(odata, x3, x5, x4, x2, x1, x0);
}

void Task1652::Task_local::compute() {
  const Index x4 = b(0);
  const Index x5 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1756
  std::unique_ptr<double[]> odata = out()->move_block(x4, x5, a2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x5, a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x5, a2, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size());
    // tensor label: I1757
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x5, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x5, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x5.size(), a2.size(), c1.size());
    dgemm_("T", "N", x3.size(), x4.size()*x5.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x4.size(), x5.size(), a2.size());
  out()->put_block(odata, x4, x5, a2, x3);
}

void Task1653::Task_local::compute() {
  const Index x4 = b(0);
  const Index x5 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1757
  std::unique_ptr<double[]> odata = out()->move_block(x4, x5, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), x5.size(), x4.size());
  }
  out()->put_block(odata, x4, x5, a2, c1);
}

void Task1654::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1417
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, x4, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, x4, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, x4, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());
    // tensor label: I2036
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<2,0,1,3,4,5,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x5, x4, x2, x1, x0);
}

void Task1655::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2036
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task1656::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma436
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x5, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x7, x0, x6, x5, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
              // tensor label: I1421
              std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, x7, x6, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, x7, x6, x5)]);
              sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), x7.size(), x6.size(), x5.size());
              dgemm_("T", "N", ci0.size(), 1, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
                     1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
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

void Task1657::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1421
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x7, x6, x5), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
    // tensor label: I1422
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, x7, x6, x5);
}

void Task1658::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1422
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task1659::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1421
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x7, x6, x5), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1784
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x6, a1, x7);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x6, a1, x7)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x6.size(), a1.size(), x7.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x5.size()*x6.size()*x7.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<2,1,0,5,4,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x5.size(), x6.size(), x7.size());
  out()->put_block(odata, x2, x1, x0, x7, x6, x5);
}

void Task1660::Task_local::compute() {
  const Index x5 = b(0);
  const Index x6 = b(1);
  const Index a1 = b(2);
  const Index x7 = b(3);
  // tensor label: I1784
  std::unique_ptr<double[]> odata = out()->move_block(x5, x6, a1, x7);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x7.size(), a1.size(), x6.size(), x5.size());
  }
  out()->put_block(odata, x5, x6, a1, x7);
}

void Task1661::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma437
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
              // tensor label: I1424
              std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, x5, x4, x3)]);
              sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), x5.size(), x4.size(), x3.size());
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

void Task1662::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: I1425
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1663::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I1425
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());
    // tensor label: I1426
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", a2.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a2);
}

void Task1664::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1426
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task1665::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1437
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, a1, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, a1, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), a1.size(), x4.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x3.size()*x5.size()*x4.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<2,1,0,4,5,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x5.size(), x4.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1666::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I1437
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
    // tensor label: I1438
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, x5, a1, x4);
}

void Task1667::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I1438
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a2.size(), x3.size());
  }
  out()->put_block(odata, a2, x3);
}

void Task1668::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I1545
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a1, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a1, x0, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), a1.size(), x0.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<5,3,4,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1669::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I1545
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, x0, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    dscal_(x2.size()*x1.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x1, a1, x0, x2);
}

void Task1670::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I1545
  std::unique_ptr<double[]> odata = out()->move_block(x1, a1, x0, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a1, x0, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a1, x0, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I1546
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", x2.size(), x1.size()*a1.size()*x0.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), a1.size(), x0.size());
  out()->put_block(odata, x1, a1, x0, x2);
}

void Task1671::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1546
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task1672::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), x1.size(), x2.size());
    // tensor label: I1787
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, x5, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, x5, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), x5.size(), a2.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<2,1,0,5,4,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1673::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index x5 = b(2);
  const Index a2 = b(3);
  // tensor label: I1787
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4, x5, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x4, x5, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x4, x5, a2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a2.size(), a1.size());
    // tensor label: I1788
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, a1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, a1, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), a1.size(), x5.size());
    dgemm_("T", "N", a2.size(), x3.size()*x4.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x4.size(), x5.size());
  out()->put_block(odata, x3, x4, x5, a2);
}

void Task1674::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1788
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4, a1, x5);
}

void Task1675::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I1799
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x0, a1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x0, a1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x0.size(), a1.size(), x1.size());
    dgemm_("T", "N", x3.size()*x4.size()*x5.size(), x2.size()*x0.size()*x1.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size()*x4.size()*x5.size());
  }
  sort_indices<3,5,4,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1676::Task_local::compute() {
  const Index x2 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1799
  std::unique_ptr<double[]> odata = out()->move_block(x2, x0, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x0, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x0, a1, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
    // tensor label: I1800
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size());
    dgemm_("T", "N", x0.size()*a1.size()*x1.size(), x2.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size()*a1.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size(), x1.size(), x2.size());
  out()->put_block(odata, x2, x0, a1, x1);
}

void Task1677::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  // tensor label: I1800
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x2);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a2.size(), x2.size());
  }
  out()->put_block(odata, a2, x2);
}

void Task1678::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1907
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, x5, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, x5, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x5.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<2,1,0,4,3,5,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1679::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I1907
  std::unique_ptr<double[]> odata = out()->move_block(x4, a1, x5, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    dscal_(x3.size()*x4.size()*a1.size()*x5.size(), e0_, i0data.get(), 1);
    sort_indices<2,1,0,3,1,1,-1,4>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x4, a1, x5, x3);
}

void Task1680::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I1907
  std::unique_ptr<double[]> odata = out()->move_block(x4, a1, x5, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, a1, x5, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x5, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size());
    // tensor label: I1908
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x4, a1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x4, a1, x5)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x4.size(), a1.size(), x5.size());
    dgemm_("T", "N", x3.size(), x4.size()*a1.size()*x5.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x4.size(), a1.size(), x5.size());
  out()->put_block(odata, x4, a1, x5, x3);
}

void Task1681::Task_local::compute() {
  const Index a2 = b(0);
  const Index x4 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1908
  std::unique_ptr<double[]> odata = out()->move_block(a2, x4, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), a2.size());
  }
  out()->put_block(odata, a2, x4, a1, x5);
}

void Task1682::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I2033
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1683::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2033
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task1684::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1424
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I2087
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, a1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, a1, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), a1.size(), x5.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<2,1,0,5,4,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
  out()->put_block(odata, x2, x1, x0, x5, x4, x3);
}

void Task1685::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I2087
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4, a1, x5);
}

void Task1686::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma438
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x0, x2, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          // tensor label: I1428
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
          sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
          dgemm_("T", "N", ci0.size(), 1, x3.size()*x2.size()*x1.size()*x0.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x3.size()*x2.size()*x1.size()*x0.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1687::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1429
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1688::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I1429
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I1430
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task1689::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1430
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task1690::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I1429
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I1434
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task1691::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1434
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task1692::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
    // tensor label: I1495
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1693::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  // tensor label: I1495
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
      // tensor label: I1496
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", 1, a3.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size());
  out()->put_block(odata, a3, x0);
}

void Task1694::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1496
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1695::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
    // tensor label: I1499
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1696::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1499
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
      // tensor label: I1500
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", 1, a1.size()*x0.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size());
  out()->put_block(odata, a1, x0);
}

void Task1697::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1500
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1698::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      // tensor label: I1541
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a3.size()*a1.size(),
             1.0, i0data_sorted, a3.size()*a1.size(), i1data_sorted, a3.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1699::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1541
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a1, x0, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());
    // tensor label: I1542
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*a1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, a1, x0, x1);
}

