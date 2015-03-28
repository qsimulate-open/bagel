//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks15.cc
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

#include <src/smith/CASPT2_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task700::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I691
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-4,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task701::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  // tensor label: I690
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I694
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", c4.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size());
  out()->put_block(odata, x0, c4);
}

void Task702::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I694
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task703::Task_local::compute() {
  const Index c4 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
  {
    // tensor label: I701
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c4, c2);
}

void Task704::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I701
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
        // tensor label: I702
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,0,2,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c2.size(), a3.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*a1.size()*x1.size(), i1data_sorted, a3.size()*a1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c2.size());
  out()->put_block(odata, c2, c4);
}

void Task705::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I702
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I703
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task706::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I703
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task707::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I701
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I705
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c2.size(), a3.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*a1.size()*x1.size(), i1data_sorted, a3.size()*a1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c2.size());
  out()->put_block(odata, c2, c4);
}

void Task708::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I705
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I706
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task709::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I706
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task710::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
  {
    // tensor label: I707
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a1.size(), a4.size());
  }
  out()->put_block(odata, a1, a4);
}

void Task711::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: I707
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
        // tensor label: I708
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a1.size(), a3.size()*c2.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*x1.size(), i1data_sorted, a3.size()*c2.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a1.size());
  out()->put_block(odata, a1, a4);
}

void Task712::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I708
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I709
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task713::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I709
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task714::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: I707
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
        // tensor label: I711
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a1.size(), a3.size()*c2.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*x1.size(), i1data_sorted, a3.size()*c2.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a1.size());
  out()->put_block(odata, a1, a4);
}

void Task715::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I711
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I712
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task716::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I712
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task717::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4);
  {
    // tensor label: I713
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a3.size(), a4.size());
  }
  out()->put_block(odata, a3, a4);
}

void Task718::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  // tensor label: I713
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
        // tensor label: I714
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a3.size(), c2.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x1.size(), i1data_sorted, c2.size()*a1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a3.size());
  out()->put_block(odata, a3, a4);
}

void Task719::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I714
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I715
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task720::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I715
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task721::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  // tensor label: I713
  std::unique_ptr<double[]> odata = out()->move_block(a3, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
        // tensor label: I717
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a3.size(), c2.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x1.size(), i1data_sorted, c2.size()*a1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a3.size());
  out()->put_block(odata, a3, a4);
}

void Task722::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I717
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I718
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task723::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I718
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task724::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
  {
    // tensor label: I719
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c2.size());
  }
  out()->put_block(odata, x1, c2);
}

void Task725::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  // tensor label: I719
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I720
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a3, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a3, x0, x1)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), a3.size(), x0.size(), x1.size());
        dgemm_("T", "N", c2.size(), x1.size(), a1.size()*a3.size()*x0.size(),
               1.0, i0data_sorted, a1.size()*a3.size()*x0.size(), i1data_sorted, a1.size()*a3.size()*x0.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size());
  out()->put_block(odata, x1, c2);
}

void Task726::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I720
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a3, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a3, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I721
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a3.size());
  out()->put_block(odata, a1, a3, x0, x1);
}

void Task727::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  // tensor label: I721
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, a3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a3.size());
  }
  out()->put_block(odata, x3, a1, x2, a3);
}

void Task729::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0);
  {
    // tensor label: I734
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c1.size(), x0.size());
  }
  out()->put_block(odata, c1, x0);
}

void Task730::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I734
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma12
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I735
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->put_block(odata, c1, x0);
}

void Task731::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I735
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task732::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  {
    // tensor label: I736
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size());
  }
  out()->put_block(odata, c1, a2);
}

void Task733::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I736
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma38
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I737
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

void Task734::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I737
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, c1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), a2.size(), c1.size(), x0.size());
  }
  out()->put_block(odata, x1, a2, c1, x0);
}

void Task735::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I736
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I739
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, a2, c1);
}

void Task736::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I739
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task737::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  {
    // tensor label: I740
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), a1.size());
  }
  out()->put_block(odata, x0, a1);
}

void Task738::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: I740
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
        // tensor label: I741
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

void Task739::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I741
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task741::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, c1, x0);
  {
    // tensor label: I742
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x0, x1);
    sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, c2, x1, c1, x0);
}

void Task742::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I742
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma92
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I743
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), c2.size());
  out()->put_block(odata, c1, c2, x0, x1);
}

void Task743::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I743
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), x3.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, c1, x3, c2, x2);
}

void Task744::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x0, x1);
  {
    // tensor label: I744
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c1, x2, x0, x1);
}

void Task745::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I744
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I745
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task746::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I745
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x4.size(), c1.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c1, x3);
}

void Task747::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0, c1, a2);
  {
    // tensor label: I746
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    sort_indices<0,3,2,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x0.size());
  }
  out()->put_block(odata, c3, x0, c1, a2);
}

void Task748::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I746
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I747
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c3, a2, c1, x0);
}

void Task749::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I747
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x1);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, c3, a2, c1, x1);
}

#endif
