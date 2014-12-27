//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks17.cc
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


#include <src/smith/CAS_test_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task800::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index c1 = b(4);
  const Index c2 = b(5);
  // tensor label: I885
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0, c1, c2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      // tensor label: I886
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x5, x3, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x5, x3, x4, x1, x0)]);
      sort_indices<3,1,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*c2.size(), x2.size()*x3.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x2.size(), x3.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x3, x1, x0, c1, c2);
}

void Task801::Task_local::compute() {
  const Index x2 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I886
  std::unique_ptr<double[]> odata = out()->move_block(x2, x5, x3, x4, x1, x0);
  {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x5, x3, x4, x1, x0);
}

void Task802::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I884
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
          // tensor label: I1041
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1, a1, c2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1, a1, c2)]);
          sort_indices<0,4,3,2,1,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c2.size());
          dgemm_("T", "N", 1, x3.size()*c2.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size());
  out()->put_block(odata, x3, c2);
}

void Task803::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I1041
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I1042
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x3, x4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x3, x4, x2, x1)]);
      sort_indices<3,0,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x3.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1, a1, c2);
}

void Task804::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1042
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x3, x4, x2, x1);
  {
    // tensor label: Gamma56
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x3, x4, x2, x1);
}

void Task805::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I1041
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1, a1, c2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1045
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x0, x2, x1)]);
      sort_indices<1,0,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,2,4,5,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1, a1, c2);
}

void Task806::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1045
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma57
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}

void Task807::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4);
  {
    // tensor label: I887
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4);
}

void Task808::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I887
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& c1 : *range_[0]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
          // tensor label: I888
          std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x4, x3, x1, x0, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x4, x3, x1, x0, c1)]);
          sort_indices<4,3,5,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x4.size(), x3.size(), x1.size(), x0.size(), c1.size());
          dgemm_("T", "N", 1, x4.size()*x3.size(), x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*c1.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->put_block(odata, x4, x3);
}

void Task809::Task_local::compute() {
  const Index x2 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  const Index c1 = b(5);
  // tensor label: I888
  std::unique_ptr<double[]> odata = out()->move_block(x2, x4, x3, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x4, x3, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x4, x3, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        // tensor label: I889
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
        sort_indices<3,1,0,2,4,5,6,7,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x4.size()*x3.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,4,5,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x4, x3, x1, x0, c1);
}

void Task810::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x1 = b(6);
  const Index x0 = b(7);
  // tensor label: I889
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, x2, x5, x4, x3, x1, x0);
  {
    // tensor label: Gamma273
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,4>(i0data, odata, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x7, x6, x2, x5, x4, x3, x1, x0);
}

void Task811::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I887
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
          // tensor label: I1047
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x4, x3, x2, x1, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x4, x3, x2, x1, a1)]);
          sort_indices<0,5,4,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x4.size(), x3.size(), x2.size(), x1.size(), a1.size());
          dgemm_("T", "N", 1, x4.size()*x3.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->put_block(odata, x4, x3);
}

void Task812::Task_local::compute() {
  const Index x0 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index a1 = b(5);
  // tensor label: I1047
  std::unique_ptr<double[]> odata = out()->move_block(x0, x4, x3, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x4, x3, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x4, x3, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        // tensor label: I1048
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
        sort_indices<3,2,0,1,4,5,6,7,0,1,1,1>(i1data, i1data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x4.size()*x3.size()*x2.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,5,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x4, x3, x2, x1, a1);
}

void Task813::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: I1048
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x4, x3, x2, x1);
  {
    // tensor label: Gamma326
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,4>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x4, x3, x2, x1);
}

void Task814::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1);
  {
    // tensor label: I890
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), c1.size());
  }
  out()->put_block(odata, c2, c1);
}

void Task815::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  // tensor label: I890
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        // tensor label: I891
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c2)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c2.size());
        dgemm_("T", "N", c1.size(), c2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size());
  out()->put_block(odata, c2, c1);
}

void Task816::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  // tensor label: I891
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I892
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<3,1,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c2.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2);
}

void Task817::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I892
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task818::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  {
    // tensor label: I893
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), a3.size());
  }
  out()->put_block(odata, c2, a3);
}

void Task819::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I893
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& c1 : *range_[0]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
          // tensor label: I894
          std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c2, a3, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c2, a3, c1)]);
          sort_indices<2,1,5,0,3,4,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c2.size(), a3.size(), c1.size());
          dgemm_("T", "N", 1, c2.size()*a3.size(), x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*c1.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->put_block(odata, c2, a3);
}

void Task820::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  // tensor label: I894
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
    // tensor label: I895
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*c1.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size()*a3.size()*c1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2, a3, c1);
}

void Task821::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I895
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task822::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  // tensor label: I894
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    // tensor label: I898
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size()*a3.size()*c2.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a3.size()*c2.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), c2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2, a3, c1);
}

void Task823::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I898
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task824::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I893
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
          // tensor label: I1053
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a3, c2, a1)]);
          sort_indices<0,5,2,1,3,4,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, a3.size()*c2.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size());
  out()->put_block(odata, c2, a3);
}

void Task825::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I1053
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1054
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}

void Task826::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1054
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task827::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I1053
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1057
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}

void Task828::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1057
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task829::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2);
  {
    // tensor label: I899
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x3.size(), a2.size());
  }
  out()->put_block(odata, x3, a2);
}

void Task830::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I899
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& c1 : *range_[0]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
          // tensor label: I900
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0, a2, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0, a2, c1)]);
          sort_indices<3,2,5,1,0,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size(), a2.size(), c1.size());
          dgemm_("T", "N", 1, x3.size()*a2.size(), x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*c1.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->put_block(odata, x3, a2);
}

void Task831::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I900
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0, a2, c1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I901
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x3, x2, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x3, x2, x4, x1, x0)]);
      sort_indices<3,0,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x3.size()*x2.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0, a2, c1);
}

void Task832::Task_local::compute() {
  const Index x5 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I901
  std::unique_ptr<double[]> odata = out()->move_block(x5, x3, x2, x4, x1, x0);
  {
    // tensor label: Gamma9
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x3, x2, x4, x1, x0);
}

void Task833::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I900
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0, a2, c1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I904
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
      sort_indices<1,0,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x2.size()*x3.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,2,4,5,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x2.size(), x3.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0, a2, c1);
}

void Task834::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I904
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task835::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I899
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
          // tensor label: I1059
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1, a1, a2)]);
          sort_indices<0,4,3,2,1,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size(), a1.size(), a2.size());
          dgemm_("T", "N", 1, x3.size()*a2.size(), x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*x2.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x2.size()*x1.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->put_block(odata, x3, a2);
}

void Task836::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index a2 = b(5);
  // tensor label: I1059
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1, a1, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      // tensor label: I1060
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
      sort_indices<2,0,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x3.size()*x2.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1, a1, a2);
}

void Task837::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1060
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma59
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}

void Task838::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1);
  {
    // tensor label: I905
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1);
}

void Task839::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  // tensor label: I905
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I906
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, c3)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), c3.size());
        dgemm_("T", "N", a2.size(), x1.size(), x0.size()*c1.size()*c3.size(),
               1.0, i0data_sorted, x0.size()*c1.size()*c3.size(), i1data_sorted, x0.size()*c1.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size());
  out()->put_block(odata, x1, a2);
}

void Task840::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I906
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      // tensor label: I907
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x0, x2)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, c3);
}

void Task841::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I907
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x0, x2);
  {
    // tensor label: Gamma3
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, x1, x3, x0, x2);
}

void Task842::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  {
    // tensor label: I908
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), c1.size());
  }
  out()->put_block(odata, a2, c1);
}

void Task843::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I908
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I909
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size());
      dgemm_("T", "N", a2.size()*c1.size(), 1, x0.size()*c3.size(),
             1.0, i0data_sorted, x0.size()*c3.size(), i1data_sorted, x0.size()*c3.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task844::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: I909
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c3, x1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        // tensor label: I910
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        dgemm_("T", "N", c3.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
  out()->put_block(odata, x0, c3);
}

void Task845::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I910
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task846::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3);
  {
    // tensor label: I911
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), a2.size());
  }
  out()->put_block(odata, a2, c3);
}

void Task847::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I911
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I912
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, x0.size()*c1.size(),
             1.0, i0data_sorted, x0.size()*c1.size(), i1data_sorted, x0.size()*c1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size());
  out()->put_block(odata, c3, a2);
}

void Task848::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  // tensor label: I912
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, x1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        // tensor label: I913
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        dgemm_("T", "N", c1.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size());
  out()->put_block(odata, x0, c1);
}

void Task849::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I913
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

