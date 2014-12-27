//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks14.cc
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


#include <src/smith/CAS_test_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task650::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I643
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I644
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c2);
}

void Task651::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I644
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task652::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I643
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I652
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c2);
}

void Task653::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I652
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task654::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I647
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, c2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), c2.size());
    dgemm_("T", "N", a3.size(), x0.size()*a1.size()*c2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task655::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I647
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I648
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, c2);
}

void Task656::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I648
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task657::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I647
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I656
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x1, a1, c2);
}

void Task658::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I656
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task659::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I659
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size());
  dgemm_("T", "N", c2.size()*a1.size(), x0.size()*a3.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), a3.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task660::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  // tensor label: I659
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        // tensor label: I660
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a3.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size());
  out()->put_block(odata, x0, a3);
}

void Task661::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I660
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task662::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I663
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size());
  dgemm_("T", "N", c2.size()*a3.size(), x0.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x0.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task663::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: I663
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        // tensor label: I664
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size());
  out()->put_block(odata, x0, a1);
}

void Task664::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I664
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task665::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());
      // tensor label: I667
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
      sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());
      dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task666::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);
  // tensor label: I667
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
  // tensor label: I668
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}

void Task667::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I668
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task668::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);
  // tensor label: I667
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
  // tensor label: I672
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}

void Task669::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I672
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task670::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I675
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task671::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I675
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma81
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}

void Task672::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I678
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task673::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I678
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma81
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}

void Task674::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());
    // tensor label: I681
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());
    dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task675::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);
  // tensor label: I681
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I682
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}

void Task676::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I682
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task677::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c4 = b(2);
  const Index a1 = b(3);
  // tensor label: I681
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I686
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c4, a1);
}

void Task678::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I686
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task679::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I689
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());
    dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task680::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I689
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I690
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}

void Task681::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I690
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task682::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  // tensor label: I689
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I694
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a3);
}

void Task683::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I694
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task684::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I697
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());
    dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task685::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I697
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I698
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}

void Task686::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I698
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task687::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I697
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I702
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c2, a1);
}

void Task688::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I702
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task689::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());
    // tensor label: I705
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());
    dgemm_("T", "N", c2.size(), x0.size()*a1.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a1.size(), a3.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task690::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index a3 = b(3);
  // tensor label: I705
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      // tensor label: I706
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, a3);
}

void Task691::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I706
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task692::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I755
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task693::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I755
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task694::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I758
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task695::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I758
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task696::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I813
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task697::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I813
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task698::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I816
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task699::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I816
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

