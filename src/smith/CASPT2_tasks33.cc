//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks33.cc
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


#include <src/smith/CASPT2_tasks33.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1600::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1887
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x1, a4), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I1888
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
    dgemm_("T", "N", a4.size(), a3.size()*c2.size()*x1.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), a3.size(), c2.size(), x1.size());
  out()->put_block(odata, a3, c2, x1, a4);
}

void Task1601::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1888
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1602::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a4.size());
        // tensor label: I1891
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x1, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x1.size(), a4.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a4.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a4.size(), i1data_sorted, a3.size()*c2.size()*a4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1603::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1891
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x1, a4), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I1892
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
    dgemm_("T", "N", a4.size(), a3.size()*c2.size()*x1.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), a3.size(), c2.size(), x1.size());
  out()->put_block(odata, a3, c2, x1, a4);
}

void Task1604::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1892
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1605::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a1.size());
        // tensor label: I1895
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x1, a4)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x1.size(), a4.size());
        dgemm_("T", "N", x0.size(), x1.size(), c2.size()*a1.size()*a4.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*a4.size(), i1data_sorted, c2.size()*a1.size()*a4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1606::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1895
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, a4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I1896
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
    dgemm_("T", "N", a4.size(), c2.size()*a1.size()*x1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x1.size());
  out()->put_block(odata, c2, a1, x1, a4);
}

void Task1607::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1896
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1608::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a4.size());
        // tensor label: I1899
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x1, a4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x1.size(), a4.size());
        dgemm_("T", "N", x0.size(), x1.size(), c2.size()*a1.size()*a4.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*a4.size(), i1data_sorted, c2.size()*a1.size()*a4.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1609::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x1 = b(2);
  const Index a4 = b(3);
  // tensor label: I1899
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x1, a4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I1900
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
    dgemm_("T", "N", a4.size(), c2.size()*a1.size()*x1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x1.size());
  out()->put_block(odata, c2, a1, x1, a4);
}

void Task1610::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1900
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1611::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I1949
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1612::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1949
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    dscal_(a3.size()*c2.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1613::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I1952
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1614::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1952
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    dscal_(a3.size()*c2.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1615::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I1985
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1616::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1985
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    dscal_(a3.size()*c2.size()*a1.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1617::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I1988
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1618::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1988
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    dscal_(a3.size()*c2.size()*a1.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1619::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I2039
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1620::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2039
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1621::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I2042
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1622::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2042
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1623::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I2093
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1624::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I2093
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1625::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I2096
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->put_block(odata, x1, x0);
}

void Task1626::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I2096
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task1627::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
      // tensor label: I2105
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", 1, x1.size()*x0.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1628::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2105
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, x0);
    sort_indices<0,2,1,3,1,1,-1,1>(i1data, odata, x1.size(), a1.size(), c2.size(), x0.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task1629::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1343
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: h1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: I2108
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x1.size()*x0.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task1630::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I2108
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,3,1,0,1,1,2,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1631::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma429
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x2, x4, x3, x1, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I1393
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, x5, x4, x3)]);
              sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
              dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
                     1.0, i0data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(), i1data_sorted, x1.size()*x0.size()*x2.size()*x5.size()*x4.size()*x3.size(),
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

void Task1632::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1393
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x2, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: I1394
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->put_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task1633::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I1394
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, x2), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I1395
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
    dgemm_("T", "N", x2.size(), x1.size()*x0.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, x1, x0, a2, x2);
}

void Task1634::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1395
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task1635::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1393
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x2, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
    // tensor label: I1780
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, x0, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), x0.size(), x1.size());
    dgemm_("T", "N", x3.size()*x4.size()*x5.size(), x2.size()*x0.size()*x1.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size()*x4.size()*x5.size());
  }
  sort_indices<5,4,3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task1636::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1780
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, x0, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
    // tensor label: I1781
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size());
    dgemm_("T", "N", a1.size()*x0.size()*x1.size(), x2.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, a1.size()*x0.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, x2, a1, x0, x1);
}

void Task1637::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  // tensor label: I1781
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x2.size(), c2.size());
  }
  out()->put_block(odata, x2, c2);
}

void Task1638::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I1393
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x2, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, x2, a1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), x2.size(), a1.size());
    // tensor label: I2090
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x4, a1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x4, a1, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x4.size(), a1.size(), x5.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size()*x4.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<1,0,2,5,4,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
  out()->put_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task1639::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I2090
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4, a1, x5);
}

void Task1640::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma434
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x4, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x0, x3, x4, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
              // tensor label: I1413
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, x4, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, x4, x2, x1, x0)]);
              sort_indices<1,5,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), x4.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", ci0.size(), 1, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(),
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

void Task1641::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, x4, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, x4, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, x4, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1414
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, a1, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, a1, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), a1.size(), x4.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x3.size()*x5.size()*x4.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x5.size(), x4.size());
  out()->put_block(odata, x3, x5, x4, x2, x1, x0);
}

void Task1642::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I1414
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    // tensor label: I1415
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, x5, a1, x4);
}

void Task1643::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1415
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task1644::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1413
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, x4, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, x4, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, x4, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1702
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, x5, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, x5, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x5.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x4.size()*x5.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<5,4,3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x4.size(), x5.size(), x3.size());
  out()->put_block(odata, x3, x5, x4, x2, x1, x0);
}

void Task1645::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  // tensor label: I1702
  std::unique_ptr<double[]> odata = out()->move_block(x4, a1, x5, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, a1, x5, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x5, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size());
    // tensor label: I1703
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c2, a1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c2, a1, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c2.size(), a1.size(), x5.size());
    dgemm_("T", "N", x3.size(), x4.size()*a1.size()*x5.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x3.size(), x4.size(), a1.size(), x5.size());
  out()->put_block(odata, x4, a1, x5, x3);
}

void Task1646::Task_local::compute() {
  const Index x4 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1703
  std::unique_ptr<double[]> odata = out()->move_block(x4, c2, a1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x5.size(), a1.size(), c2.size(), x4.size());
  }
  out()->put_block(odata, x4, c2, a1, x5);
}

void Task1647::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x0 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma435
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
              // tensor label: I1417
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, x4, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, x4, x2, x1, x0)]);
              sort_indices<1,2,0,5,3,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), x4.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", ci0.size(), 1, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x3.size()*x5.size()*x4.size()*x2.size()*x1.size()*x0.size(),
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

void Task1648::Task_local::compute() {
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
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1418
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x5, x4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x5.size(), x4.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x3.size()*x5.size()*x4.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size(), x5.size(), x4.size());
  out()->put_block(odata, x3, x5, x4, x2, x1, x0);
}

void Task1649::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I1418
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
    // tensor label: I1419
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", a1.size()*x5.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, a1.size()*x5.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

