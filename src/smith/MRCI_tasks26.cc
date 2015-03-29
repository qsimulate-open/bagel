//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks26.cc
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

#include <src/smith/MRCI_tasks26.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1250::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a3.size());
      // tensor label: I1608
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1251::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1608
  std::unique_ptr<double[]> odata = out()->move_block(c2, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1609
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->put_block(odata, c2, x5, x0, x4);
}

void Task1252::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1609
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, c2, x3, x2, x1);
}

void Task1253::Task_local::compute() {
  const Index c2 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I1608
  std::unique_ptr<double[]> odata = out()->move_block(c2, x5, x0, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma526
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: I1612
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c2, x1)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c2.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c2.size());
  out()->put_block(odata, c2, x5, x0, x4);
}

void Task1254::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1612
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c2, x1);
}

void Task1255::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1614
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1256::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1614
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1615
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a3.size(), a4.size());
  out()->put_block(odata, a3, a4, x0, x1);
}

void Task1257::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1615
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, a4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), a3.size(), x2.size(), a4.size());
  }
  out()->put_block(odata, x3, a3, x2, a4);
}

void Task1258::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c2.size(), a3.size());
      // tensor label: I1617
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1259::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1617
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1618
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->put_block(odata, a1, a4, x0, x1);
}

void Task1260::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1618
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a4.size());
  }
  out()->put_block(odata, x3, a1, x2, a4);
}

void Task1261::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, a4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1, a4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a1.size());
      // tensor label: I1620
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1262::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1620
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1621
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), a4.size());
  out()->put_block(odata, a3, a4, x1, x0);
}

void Task1263::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1621
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, a4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a3.size(), x2.size(), a4.size());
  }
  out()->put_block(odata, x3, a3, x2, a4);
}

void Task1264::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, a4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1, a4, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size(), a4.size(), a3.size());
      // tensor label: I1623
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), a4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1265::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1623
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1624
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a4.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a4.size());
  out()->put_block(odata, a1, a4, x0, x1);
}

void Task1266::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a4 = b(3);
  // tensor label: I1624
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a4.size());
  }
  out()->put_block(odata, x3, a1, x2, a4);
}

void Task1267::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma562
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I1705
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1268::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1705
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c2, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x3, a3, c2, a1);
}

void Task1269::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma564
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size());
    // tensor label: I1709
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1270::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1709
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, c2, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x5, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,-1,1>(i1data, odata, x5.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x5, a3, c2, a1);
}

void Task1271::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
  {
    // tensor label: I239
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x1, a2);
    sort_indices<2,3,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x1.size(), a2.size());
  }
  {
    // tensor label: I239
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x1, x0, a1);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, a2.size(), x1.size(), x0.size(), a1.size());
  }
  out()->put_block(odata, x1, a2, x0, a1);
}

void Task1272::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I240
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a2.size(), a1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), a1.size(), x0.size(), x1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1273::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I240
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I241
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task1274::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I241
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1275::Task_local::compute() {
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
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I243
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a1, a2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1276::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
    // tensor label: I244
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", x3.size()*a1.size()*a2.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), a2.size(), x2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1277::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I244
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task1278::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
    // tensor label: I247
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", x3.size()*a1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*x2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1279::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I247
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task1280::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c4.size(), a3.size());
      // tensor label: I1654
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1281::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1654
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, x2, c4, a3, a2);
}

void Task1282::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I1657
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1283::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I1657
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a2, a4, c3);
}

void Task1284::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I1660
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*a1.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1285::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I1660
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x2.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a2, a4, c3);
}

void Task1286::Task_local::compute() {
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
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), x4.size());
        // tensor label: I1626
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x5, x0, x4, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x5, x0, x4, x1)]);
        sort_indices<2,1,4,0,3,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x5.size(), x0.size(), x4.size(), x1.size());
        dgemm_("T", "N", a1.size(), a2.size()*x0.size()*x1.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1287::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x4 = b(4);
  const Index x1 = b(5);
  // tensor label: I1626
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x5, x0, x4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x5, x0, x4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x5, x0, x4, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma531
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x2, x4, x3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x2, x4, x3, x1)]);
      sort_indices<2,4,0,1,3,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x2.size(), x4.size(), x3.size(), x1.size());
      // tensor label: I1627
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, c3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x0.size()*x4.size()*x1.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x0.size()*x4.size()*x1.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), x1.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c3, x5, x0, x4, x1);
}

void Task1288::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1627
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a2.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x3, a2, x2, c3);
}

void Task1289::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I1629
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x5, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x5, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", a1.size(), a2.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1290::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1629
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma532
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x1, x2, x0)]);
      sort_indices<2,4,0,1,3,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1630
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, c3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c3, x5, x4, x1, x0);
}

void Task1291::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1630
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, x2, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x3, a2, x2, c3);
}

void Task1292::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I1632
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x7, x0, x6, x5, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x7, x0, x6, x5, x1)]);
        sort_indices<1,3,4,0,2,5,0,1,1,1>(i1data, i1data_sorted, a2.size(), x7.size(), x0.size(), x6.size(), x5.size(), x1.size());
        dgemm_("T", "N", a1.size(), a2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1293::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1632
  std::unique_ptr<double[]> odata = out()->move_block(a2, x7, x0, x6, x5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma533
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x1, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x1, x3, x2)]);
        sort_indices<4,6,7,0,1,2,3,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size());
        // tensor label: I1633
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a2, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a2, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), a2.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->put_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task1294::Task_local::compute() {
  const Index x4 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1633
  std::unique_ptr<double[]> odata = out()->move_block(x4, a2, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x4, a2, x3, x2);
}

void Task1295::Task_local::compute() {
  const Index a2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I1632
  std::unique_ptr<double[]> odata = out()->move_block(a2, x7, x0, x6, x5, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x7, x0, x6, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x7, x0, x6, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma349
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
        sort_indices<4,5,6,0,1,2,3,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1636
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x2, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x2, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x2.size(), a2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x5.size()*x1.size(), a2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), a2.size());
  out()->put_block(odata, a2, x7, x0, x6, x5, x1);
}

void Task1296::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I1636
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, x2, a2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), x3.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, x4, x3, x2, a2);
}

void Task1297::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, a3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), a3.size(), a2.size());
      // tensor label: I1638
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1, x2, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, a1, x0, x1, a2);
}

void Task1298::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1638
  std::unique_ptr<double[]> odata = out()->move_block(a3, x1, x2, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma471
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x4, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x4, x3, x2, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: I1639
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x1, x2, x0);
}

void Task1299::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1639
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

#endif
