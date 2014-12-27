//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks20.cc
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


#include <src/smith/CASPT2_tasks20.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task950::Task_local::compute() {
  const Index a2 = b(0);
  const Index a3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a3);
  {
    // tensor label: I1013
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a2, a3);
}

void Task951::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I1013
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
        // tensor label: I1014
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a3, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a3, c1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a3.size(), c1.size());
        dgemm_("T", "N", a2.size(), a3.size(), x1.size()*x0.size()*c1.size(),
               1.0, i0data_sorted, x1.size()*x0.size()*c1.size(), i1data_sorted, x1.size()*x0.size()*c1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a3.size());
  out()->put_block(odata, a3, a2);
}

void Task952::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I1014
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
      // tensor label: I1015
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c1);
}

void Task953::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1015
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task954::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c1 = b(3);
  // tensor label: I1014
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I1024
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a3.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c1);
}

void Task955::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1024
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task956::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I1013
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
        // tensor label: I1160
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
        sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());
        dgemm_("T", "N", a2.size(), a3.size(), x0.size()*x1.size()*a1.size(),
               1.0, i0data_sorted, x0.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x1.size()*a1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a3.size());
  out()->put_block(odata, a3, a2);
}

void Task957::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);
  // tensor label: I1160
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      // tensor label: I1161
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a3);
}

void Task958::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1161
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task959::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1);
  {
    // tensor label: I1025
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x2.size(), c1.size());
  }
  out()->put_block(odata, x2, c1);
}

void Task960::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  // tensor label: I1025
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
        // tensor label: I1026
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a2)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
        dgemm_("T", "N", c1.size(), x2.size(), x1.size()*x0.size()*a2.size(),
               1.0, i0data_sorted, x1.size()*x0.size()*a2.size(), i1data_sorted, x1.size()*x0.size()*a2.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size());
  out()->put_block(odata, x2, c1);
}

void Task961::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I1026
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I1027
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<3,2,0,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a2);
}

void Task962::Task_local::compute() {
  const Index x5 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1027
  std::unique_ptr<double[]> odata = out()->move_block(x5, x2, x4, x3, x1, x0);
  {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x4, x3, x1, x0);
}

void Task963::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2);
  {
    // tensor label: I1049
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a2.size(), a1.size());
  }
  out()->put_block(odata, a1, a2);
}

void Task964::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  // tensor label: I1049
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        // tensor label: I1050
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a2)]);
        sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a2.size());
        dgemm_("T", "N", a1.size(), a2.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size());
  out()->put_block(odata, a2, a1);
}

void Task965::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1050
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        // tensor label: I1051
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<3,2,0,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a2.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a2);
}

void Task966::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1051
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma59
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}

void Task967::Task_local::compute() {
  const Index a2 = b(0);
  const Index x0 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, x0);
  {
    // tensor label: I1061
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x0.size(), a2.size());
  }
  out()->put_block(odata, a2, x0);
}

void Task968::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  // tensor label: I1061
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1062
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1, a4, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1, a4, c3)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size(), a4.size(), c3.size());
        dgemm_("T", "N", a2.size(), x0.size(), c1.size()*a4.size()*c3.size(),
               1.0, i0data_sorted, c1.size()*a4.size()*c3.size(), i1data_sorted, c1.size()*a4.size()*c3.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size());
  out()->put_block(odata, x0, a2);
}

void Task969::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I1062
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I1063
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c1, a4, c3);
}

void Task970::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1063
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task971::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a4, x0);
  {
    // tensor label: I1064
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x0.size(), a4.size());
  }
  out()->put_block(odata, a4, x0);
}

void Task972::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I1064
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1065
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1, a2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1, a2, c3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size(), a2.size(), c3.size());
        dgemm_("T", "N", a4.size(), x0.size(), c1.size()*a2.size()*c3.size(),
               1.0, i0data_sorted, c1.size()*a2.size()*c3.size(), i1data_sorted, c1.size()*a2.size()*c3.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->put_block(odata, x0, a4);
}

void Task973::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1065
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I1066
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c1, a2, c3);
}

void Task974::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I1066
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task975::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: I1070
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task976::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I1070
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I1071
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", a4.size()*c3.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, a4.size()*c3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c3.size());
  out()->put_block(odata, a4, c3);
}

void Task977::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I1071
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      // tensor label: I1072
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task978::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1072
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task979::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I1071
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I1078
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, a2, c1);
}

void Task980::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1078
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task981::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: I1079
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task982::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1079
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& a4 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
          // tensor label: I1080
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a4, c3, a2)]);
          sort_indices<2,5,4,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a4.size(), c3.size(), a2.size());
          dgemm_("T", "N", 1, x1.size()*x0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task983::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I1080
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1081
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task984::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1081
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task985::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I1080
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1084
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a4, c3, a2);
}

void Task986::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1084
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task987::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  {
    // tensor label: I1085
    std::unique_ptr<double[]> i0data = in(0)->get_block(c5, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c5.size(), c3.size());
  }
  out()->put_block(odata, c5, c3);
}

void Task988::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: I1085
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1086
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
        dgemm_("T", "N", c3.size(), c5.size(), c1.size()*a4.size()*a2.size(),
               1.0, i0data_sorted, c1.size()*a4.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*a2.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c5, c3);
}

void Task989::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);
  // tensor label: I1086
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
    sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}

void Task990::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
  {
    // tensor label: I1089
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a4, a5);
}

void Task991::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I1089
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1090
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
        dgemm_("T", "N", a4.size(), a5.size(), c1.size()*c3.size()*a2.size(),
               1.0, i0data_sorted, c1.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*c3.size()*a2.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a5.size());
  out()->put_block(odata, a5, a4);
}

void Task992::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1090
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
    sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
  }
  out()->put_block(odata, c1, a5, c3, a2);
}

void Task993::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  {
    // tensor label: I1093
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), c3.size());
  }
  out()->put_block(odata, x0, c3);
}

void Task994::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: I1093
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I1094
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());
        dgemm_("T", "N", c3.size(), x0.size(), a4.size()*c1.size()*a2.size(),
               1.0, i0data_sorted, a4.size()*c1.size()*a2.size(), i1data_sorted, a4.size()*c1.size()*a2.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
  out()->put_block(odata, x0, c3);
}

void Task995::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1094
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I1095
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}

void Task996::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1095
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task997::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1094
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I1098
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}

void Task998::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1098
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task999::Task_local::compute() {
  const Index a1 = b(0);
  const Index x1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a1, x1);
  {
    // tensor label: I1099
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a1.size());
  }
  out()->put_block(odata, a1, x1);
}

