//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks10.cc
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

#include <src/smith/CASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task450::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I556
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma37
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: I557
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a1.size());
  out()->put_block(odata, a1, x0, x1, x2);
}

void Task451::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I557
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task452::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  {
    // tensor label: I453
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x2.size(), a3.size());
  }
  out()->put_block(odata, x2, a3);
}

void Task453::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: I453
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
        // tensor label: I454
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c1, x3, x2)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), x2.size(), c2.size()*c1.size()*x3.size(),
               1.0, i0data_sorted, c2.size()*c1.size()*x3.size(), i1data_sorted, c2.size()*c1.size()*x3.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size());
  out()->put_block(odata, x2, a3);
}

void Task454::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I454
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma3
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I455
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, x0, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, x0, c1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), x0.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), c1.size());
  out()->put_block(odata, c2, c1, x3, x2);
}

void Task455::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I455
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, x0, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, c1.size(), x0.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, x0, c1);
}

void Task456::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: I453
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I565
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), x2.size(), c2.size()*a1.size()*x3.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x3.size(), i1data_sorted, c2.size()*a1.size()*x3.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size());
  out()->put_block(odata, x2, a3);
}

void Task457::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I566
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task458::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I566
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task459::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: I453
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I568
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), x2.size(), c2.size()*a1.size()*x3.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x3.size(), i1data_sorted, c2.size()*a1.size()*x3.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size());
  out()->put_block(odata, x2, a3);
}

void Task460::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I568
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I569
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task461::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I569
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task462::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: I453
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
        // tensor label: I607
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x3, x2)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), x2.size(), a2.size()*c1.size()*x3.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x3.size(), i1data_sorted, a2.size()*c1.size()*x3.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size());
  out()->put_block(odata, x2, a3);
}

void Task463::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I607
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I608
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, x3, x2);
}

void Task464::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I608
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task465::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  // tensor label: I453
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
        // tensor label: I610
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x3, x2)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), x2.size(), a2.size()*c1.size()*x3.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x3.size(), i1data_sorted, a2.size()*c1.size()*x3.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size());
  out()->put_block(odata, x2, a3);
}

void Task466::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I610
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I611
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x3.size()*x2.size(), a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, x3, x2);
}

void Task467::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I611
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task468::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3);
  {
    // tensor label: I456
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, c2, x3);
}

void Task469::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I456
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
        // tensor label: I457
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x3, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x3, x4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x3.size(), x4.size());
        dgemm_("T", "N", c2.size(), x3.size(), c1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c1.size()*x5.size()*x4.size(), i1data_sorted, c1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size());
  out()->put_block(odata, x3, c2);
}

void Task470::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  // tensor label: I457
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x3, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x3, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x3, x4), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x3, x4, x1, x0)]);
        sort_indices<0,4,5,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I458
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
        dgemm_("T", "N", x5.size()*x3.size()*x4.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x3.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x3.size(), x4.size(), c1.size());
  out()->put_block(odata, c1, x5, x3, x4);
}

void Task471::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I458
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task472::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I456
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
        // tensor label: I613
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x3, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x3, x4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x3.size(), x4.size());
        dgemm_("T", "N", c2.size(), x3.size(), a1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, a1.size()*x5.size()*x4.size(), i1data_sorted, a1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size());
  out()->put_block(odata, x3, c2);
}

void Task473::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  // tensor label: I613
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x3, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x3, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x3, x4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma56
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I614
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
        dgemm_("T", "N", x5.size()*x3.size()*x4.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x3.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x3.size(), x4.size(), a1.size());
  out()->put_block(odata, a1, x5, x3, x4);
}

void Task474::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I614
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task475::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I456
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I616
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", c2.size(), x3.size(), a1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, a1.size()*x5.size()*x4.size(), i1data_sorted, a1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size());
  out()->put_block(odata, x3, c2);
}

void Task476::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I616
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma57
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I617
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task477::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I617
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task478::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x3, x4);
  {
    // tensor label: I459
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x4.size(), x3.size());
  }
  out()->put_block(odata, x3, x4);
}

void Task479::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I459
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma165
              std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
              sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I460
              std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, x7, x6, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, x7, x6, x5)]);
              sort_indices<3,4,0,5,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), x7.size(), x6.size(), x5.size());
              dgemm_("T", "N", x4.size()*x3.size(), 1, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
                     1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
                     1.0, odata_sorted, x4.size()*x3.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->put_block(odata, x4, x3);
}

void Task480::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I460
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x7, x6, x5), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
    // tensor label: I461
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
    dgemm_("T", "N", x7.size()*x6.size()*x5.size(), x2.size()*x1.size()*x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, x7, x6, x5);
}

void Task481::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I461
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task482::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I459
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma218
              std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
              sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
              // tensor label: I619
              std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, x7, x6, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, x7, x6, x5)]);
              sort_indices<3,2,4,5,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), x7.size(), x6.size(), x5.size());
              dgemm_("T", "N", x4.size()*x3.size(), 1, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
                     1.0, i0data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x2.size()*x1.size()*x0.size()*x7.size()*x6.size()*x5.size(),
                     1.0, odata_sorted, x4.size()*x3.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->put_block(odata, x4, x3);
}

void Task483::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I619
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, x7, x6, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, x7, x6, x5), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
    // tensor label: I620
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

void Task484::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I620
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task485::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, c1);
  {
    // tensor label: I462
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c1.size(), c2.size());
  }
  out()->put_block(odata, c2, c1);
}

void Task486::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  // tensor label: I462
  std::unique_ptr<double[]> odata = out()->move_block(c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I463
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", c2.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), c1.size());
  out()->put_block(odata, c1, c2);
}

void Task487::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I463
  std::unique_ptr<double[]> odata = out()->move_block(c1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I464
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->put_block(odata, c1, x5, x4, x3);
}

void Task488::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I464
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task489::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  {
    // tensor label: I465
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), a3.size());
  }
  out()->put_block(odata, c2, a3);
}

void Task490::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I465
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
      // tensor label: I466
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
      dgemm_("T", "N", c2.size()*a3.size(), 1, c1.size()*x3.size(),
             1.0, i0data_sorted, c1.size()*x3.size(), i1data_sorted, c1.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->put_block(odata, c2, a3);
}

void Task491::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I466
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma7
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I467
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->put_block(odata, c1, x3);
}

void Task492::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I467
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task493::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I465
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
      // tensor label: I469
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
      dgemm_("T", "N", a3.size()*c2.size(), 1, c1.size()*x3.size(),
             1.0, i0data_sorted, c1.size()*x3.size(), i1data_sorted, c1.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size());
  out()->put_block(odata, c2, a3);
}

void Task494::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I469
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma7
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I470
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->put_block(odata, c1, x3);
}

void Task495::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I470
  std::unique_ptr<double[]> odata = out()->move_block(x2, c1, x1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), x1.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x2, c1, x1, x0);
}

void Task496::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I465
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I625
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size());
      dgemm_("T", "N", a3.size()*c2.size(), 1, a1.size()*x3.size(),
             1.0, i0data_sorted, a1.size()*x3.size(), i1data_sorted, a1.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size());
  out()->put_block(odata, c2, a3);
}

void Task497::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  // tensor label: I625
  std::unique_ptr<double[]> odata = out()->move_block(a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I626
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
        dgemm_("T", "N", x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, a1, x3);
}

void Task498::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I626
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task499::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I465
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I628
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size());
      dgemm_("T", "N", c2.size()*a3.size(), 1, a1.size()*x3.size(),
             1.0, i0data_sorted, a1.size()*x3.size(), i1data_sorted, a1.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->put_block(odata, c2, a3);
}

#endif
