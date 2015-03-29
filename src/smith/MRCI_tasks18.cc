//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks18.cc
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

#include <src/smith/MRCI_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task850::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1044
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x6, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma346
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,3,4,0,1,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1048
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task851::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1048
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, x5, x4, x3, c2);
}

void Task852::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma349
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
            // tensor label: I1056
            std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x4, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x4, x3, x7, x6, x5)]);
            sort_indices<3,4,5,1,2,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task853::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1056
  std::unique_ptr<double[]> odata = out()->move_block(a1, x4, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I1057
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, a1, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, a1, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), a1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task854::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1057
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a2.size(), a1.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a2, a1);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x4.size(), x3.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, a2, a1, x4, x3);
}

void Task855::Task_local::compute() {
  const Index a1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1056
  std::unique_ptr<double[]> odata = out()->move_block(a1, x4, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: I1105
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x5, x4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x5, x4, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x5.size(), x4.size(), x3.size());
    dgemm_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<1,4,5,0,2,3,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x4, x3, x7, x6, x5);
}

void Task856::Task_local::compute() {
  const Index a2 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1105
  std::unique_ptr<double[]> odata = out()->move_block(a2, x5, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x5, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a2.size(), x5.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, a2, x5, x4, x3);
}

void Task857::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma350
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4, x6, x5, x3, x0, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x4, x6, x5, x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I1059
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a1, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a1, x7, x6, x5)]);
            sort_indices<3,0,4,5,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task858::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1059
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a1, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, a1, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I1060
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x4, x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x4, x3, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x4.size(), x3.size(), a1.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*x3.size()*a1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, x4, x3, a1, x7, x6, x5);
}

void Task859::Task_local::compute() {
  const Index a2 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index a1 = b(3);
  // tensor label: I1060
  std::unique_ptr<double[]> odata = out()->move_block(a2, x4, x3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x4, x3, a1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a2.size(), x4.size(), x3.size(), a1.size());
  }
  out()->put_block(odata, a2, x4, x3, a1);
}

void Task860::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma351
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x5, x4, x0, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x3, x6, x5, x4, x0, x2, x1)]);
            sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
            // tensor label: I1062
            std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, x3, x7, x6, x5);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), x3.size(), x7.size(), x6.size(), x5.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task861::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I1062
  std::unique_ptr<double[]> odata = out()->move_block(x4, a1, x3, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, a1, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, a1, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: I1063
    std::unique_ptr<double[]> i1data = in(1)->get_block(x4, a1, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, a1, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), a1.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x4.size()*a1.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x4.size(), a1.size(), x3.size());
  out()->put_block(odata, x4, a1, x3, x7, x6, x5);
}

void Task862::Task_local::compute() {
  const Index x4 = b(0);
  const Index a1 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1063
  std::unique_ptr<double[]> odata = out()->move_block(x4, a1, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a1, a2, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), a1.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x4, a1, a2, x3);
}

void Task863::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma359
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x4, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x4, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1086
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, x5, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, x5, a1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x5.size(), a1.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task864::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x5 = b(2);
  const Index a1 = b(3);
  // tensor label: I1086
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, x5, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, x5, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, x5, a1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), a1.size());
      // tensor label: I1087
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, x5, a1);
}

void Task865::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1087
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, a2, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x4.size(), c3.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x4, c3, a2, x3);
}

void Task866::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma366
            std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x3, x5, x4, x2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x3, x5, x4, x2, x1)]);
            sort_indices<0,2,3,4,5,1,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x3.size(), x5.size(), x4.size(), x2.size(), x1.size());
            // tensor label: I1107
            std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x7, a1, x6);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x7, a1, x6)]);
            sort_indices<3,5,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x7.size(), a1.size(), x6.size());
            dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x0.size()*x2.size()*x1.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task867::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index a1 = b(4);
  const Index x6 = b(5);
  // tensor label: I1107
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x7, a1, x6);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x5, x4, x3, x7, a1, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, x7, a1, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
    // tensor label: I1108
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    dgemm_("T", "N", x7.size()*a1.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*a1.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), a1.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x5, x4, x3, x7, a1, x6);
}

void Task868::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1108
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, a2, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x4.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, a2, x3);
}

void Task869::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma556
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        // tensor label: I1693
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task870::Task_local::compute() {
  const Index x7 = b(0);
  const Index a1 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  // tensor label: I1693
  std::unique_ptr<double[]> odata = out()->move_block(x7, a1, x6, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x7.size(), a1.size(), x6.size(), x5.size());
  }
  out()->put_block(odata, x7, a1, x6, x5);
}

void Task871::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma557
        std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0, x8, x7, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x9, x0, x8, x7, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x0.size(), x8.size(), x7.size(), x2.size(), x1.size());
        // tensor label: I1695
        std::unique_ptr<double[]> i1data = in(1)->get_block(x9, a1, x8, x7);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x9, a1, x8, x7)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x9.size(), a1.size(), x8.size(), x7.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task872::Task_local::compute() {
  const Index x9 = b(0);
  const Index a1 = b(1);
  const Index x8 = b(2);
  const Index x7 = b(3);
  // tensor label: I1695
  std::unique_ptr<double[]> odata = out()->move_block(x9, a1, x8, x7);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, a1, x8, x7);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x9.size(), a1.size(), x8.size(), x7.size());
  }
  out()->put_block(odata, x9, a1, x8, x7);
}

void Task873::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I162
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1, a4, c3);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), a4.size(), c3.size());
  }
  {
    // tensor label: I162
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, a2, c1);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a4.size(), c3.size(), a2.size(), c1.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task874::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I163
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task875::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I163
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I164
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size());
    dgemm_("T", "N", x1.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->put_block(odata, a2, x1);
}

void Task876::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  // tensor label: I164
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x0.size(), a2.size());
  }
  out()->put_block(odata, x0, a2);
}

void Task877::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I166
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task878::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I166
  std::unique_ptr<double[]> odata = out()->move_block(a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I167
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->put_block(odata, a4, x1);
}

void Task879::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I167
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x0.size(), a4.size());
  }
  out()->put_block(odata, x0, a4);
}

void Task880::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I169
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
  dgemm_("T", "N", c3.size()*a2.size(), a4.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task881::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I169
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I170
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task882::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I170
  std::unique_ptr<double[]> odata = out()->move_block(x1, a4, c1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a4.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x1, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, a4, c1, x0);
}

void Task883::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
  // tensor label: I172
  std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
  dgemm_("T", "N", c3.size()*a4.size(), a2.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a4.size());
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task884::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I172
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I173
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task885::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I173
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, c1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x1.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,1,0,3,1,1,4,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, a2, c1, x0);
}

void Task886::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
    // tensor label: I181
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    dgemm_("T", "N", c1.size()*a4.size()*a2.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a4.size()*a2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task887::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,4,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5);
}

void Task888::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1243
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task889::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1243
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,2,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5, x1, x0);
}

void Task890::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: I1255
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task891::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c5 = b(3);
  // tensor label: I1255
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1, x0, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, x0, c5);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), x1.size(), x0.size(), c5.size());
  }
  out()->put_block(odata, c3, x1, x0, c5);
}

void Task892::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, a4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
    // tensor label: I183
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size());
    dgemm_("T", "N", c1.size()*a2.size()*a4.size(), c3.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*a4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), a4.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task893::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<0,1,1,1,-8,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5);
}

void Task894::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1246
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, c5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, c5, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task895::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1246
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5, x1, x0);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c3.size(), c5.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, c5, c3, x0);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, x1.size(), c5.size(), c3.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, c3, c5);
    sort_indices<2,3,0,1,1,1,-4,1>(i2data, odata, x1.size(), x0.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, c3, c5, x1, x0);
}

void Task896::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I183
  std::unique_ptr<double[]> odata = out()->move_block(c3, c5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: I1258
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, c5);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, c5)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), c5.size());
      dgemm_("T", "N", 1, c3.size()*c5.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c3, c5);
}

void Task897::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c5 = b(3);
  // tensor label: I1258
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1, x0, c5);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, x0, c5);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c3.size(), x1.size(), x0.size(), c5.size());
  }
  out()->put_block(odata, c3, x1, x0, c5);
}

void Task898::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
    // tensor label: I185
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task899::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4);
}

#endif
