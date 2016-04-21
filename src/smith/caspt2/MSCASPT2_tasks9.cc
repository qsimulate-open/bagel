//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_tasks9.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/caspt2/MSCASPT2_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task400::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I804
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, c3, x0)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c2.size()*a1.size()*x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), x0.size());
  out()->add_block(odata, c2, a1, x0, x1);
}

void Task401::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I371
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I854
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
      dgemm_("T", "N", x2.size()*x3.size(), x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task402::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I854
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
    // tensor label: I855
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a2.size()*c1.size()*x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c1, x0, x1);
}

void Task403::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I855
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->add_block(odata, c3, a2, c1, x0);
}

void Task404::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma124
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x3)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x3.size());
      // tensor label: I379
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x3.size(),
             1.0, i0data_sorted, x0.size()*x3.size(), i1data_sorted, x0.size()*x3.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task405::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I379
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->add_block(odata, x0, x3);
}

void Task406::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I379
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->add_block(odata, x0, x3);
}

void Task407::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I379
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
        sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x0, x3);
}

void Task408::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I379
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x3)]);
        sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x0, x3);
}

void Task409::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma126
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size());
      // tensor label: I385
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", ci0.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task410::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: l2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I386
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x1, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x1, c1)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), c1.size());
        dgemm_("T", "N", x0.size(), x1.size(), a2.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a2.size()*c3.size()*c1.size(), i1data_sorted, a2.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task411::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  // tensor label: I386
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I387
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, c3, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    dgemm_("T", "N", c1.size(), a2.size()*c3.size()*x1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x1.size());
  out()->add_block(odata, a2, c3, x1, c1);
}

void Task412::Task_local::compute() {
  const Index c4 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c4, a2, c3, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c4.size(), a2.size(), c3.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, a2, c4, x1);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), a2.size(), c4.size(), x1.size());
  }
  out()->add_block(odata, c4, a2, c3, x1);
}

void Task413::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  // tensor label: I386
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I395
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, c1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    dgemm_("T", "N", c3.size(), a2.size()*c1.size()*x1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x1.size());
  out()->add_block(odata, a2, c3, x1, c1);
}

void Task414::Task_local::compute() {
  const Index c4 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I395
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, a2, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c4, a2, c1, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c4.size(), a2.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c4, x1);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c4.size(), x1.size());
  }
  out()->add_block(odata, c4, a2, c1, x1);
}

void Task415::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  // tensor label: I386
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I399
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, c1, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    dgemm_("T", "N", a2.size(), c3.size()*c1.size()*x1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x1.size());
  out()->add_block(odata, a2, c3, x1, c1);
}

void Task416::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I399
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, c1, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, x1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a4.size(), c1.size(), x1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, c3, x1);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a4.size(), c3.size(), x1.size());
  }
  out()->add_block(odata, c3, a4, c1, x1);
}

void Task417::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
    // tensor label: I426
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task418::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I426
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->add_block(odata, x0, a4);
}

void Task419::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I426
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata.get(), out()->get_size(x0, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->add_block(odata, x0, a4);
}

void Task420::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size());
    // tensor label: I569
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task421::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I569
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,3,2,1,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), a2.size(), a4.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*c3.size()*c1.size(), i1data_sorted, a4.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->add_block(odata, a2, x1);
}

void Task422::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size());
    // tensor label: I573
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), x1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task423::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I573
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->add_block(odata, a4, x1);
}

void Task424::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: l2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I748
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x0, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x0, c1)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x0.size(), c1.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a2.size()*c3.size()*c1.size(), i1data_sorted, a2.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task425::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I748
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I749
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", c1.size(), a2.size()*c3.size()*x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task426::Task_local::compute() {
  const Index c4 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I749
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, a2, c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, a2, c3, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c4.size(), a2.size(), c3.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c3, a2, c4, x0);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c3.size(), a2.size(), c4.size(), x0.size());
  }
  out()->add_block(odata, c4, a2, c3, x0);
}

void Task427::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I748
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I757
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, c1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), c1.size(), x0.size());
    dgemm_("T", "N", c3.size(), a2.size()*c1.size()*x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task428::Task_local::compute() {
  const Index c4 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I757
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c4, a2, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c4.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c4, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c4.size(), x0.size());
  }
  out()->add_block(odata, c4, a2, c1, x0);
}

void Task429::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I748
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x0, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I761
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, c1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), c1.size(), x0.size());
    dgemm_("T", "N", a2.size(), c3.size()*c1.size()*x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->add_block(odata, a2, c3, x0, c1);
}

void Task430::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I761
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a4.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, c3, x0);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a4.size(), c3.size(), x0.size());
  }
  out()->add_block(odata, c3, a4, c1, x0);
}

void Task431::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x0.size());
    // tensor label: I788
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size());
    dgemm_("T", "N", x0.size(), x1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task432::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: I788
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", a4.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x1.size());
  out()->add_block(odata, x1, a4);
}

void Task433::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: I788
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", a4.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x1.size());
  out()->add_block(odata, x1, a4);
}

void Task434::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I931
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task435::Task_local::compute() {
  const Index a2 = b(0);
  const Index x0 = b(1);
  // tensor label: I931
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,3,2,1,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x0.size(), a2.size(), a4.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*c3.size()*c1.size(), i1data_sorted, a4.size()*c3.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a2.size());
  out()->add_block(odata, a2, x0);
}

void Task436::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size());
    // tensor label: I935
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task437::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I935
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x0.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task438::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        dscal_(x0.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task439::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
        dscal_(x0.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task440::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        dscal_(x1.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task441::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: l2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x1);
        dscal_(x1.size()*c3.size()*a2.size()*c1.size(), e0_, i1data.get(), 1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,-2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), x1.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size());
  out()->add_block(odata, x1, x0);
}

void Task442::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma132
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x1, x0, x2)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
          // tensor label: I409
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x3, x2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x3, x2)]);
          sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x3.size(), x2.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x3.size()*x2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*x3.size()*x2.size(), i1data_sorted, x0.size()*x1.size()*x3.size()*x2.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task443::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I409
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I410
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x1, x3, x2);
}

void Task444::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I410
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, x0, c3, a2, x1);
}

void Task445::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I409
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: l2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I800
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), c2.size(), x0.size(), x1.size());
      dgemm_("T", "N", x2.size()*x3.size(), x0.size()*x1.size(), a1.size()*c2.size(),
             1.0, i0data_sorted, a1.size()*c2.size(), i1data_sorted, a1.size()*c2.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->add_block(odata, x0, x1, x3, x2);
}

void Task446::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I800
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a1, c2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a1, c2, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a1.size(), c2.size(), x0.size());
    dgemm_("T", "N", x1.size(), a1.size()*c2.size()*x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), x0.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task447::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I324
  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma138
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x3, x2, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x5, x4, x1, x3, x2, x0)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
              // tensor label: I433
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, x5, x4, x3);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, x5, x4, x3)]);
              sort_indices<3,4,0,5,2,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), x5.size(), x4.size(), x3.size());
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
  out()->add_block(odata, ci0);
}

void Task448::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  // tensor label: I433
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    // tensor label: I434
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, x0, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), x0.size(), x2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), x1.size()*x0.size()*x2.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), x1.size(), x0.size(), x2.size());
  out()->add_block(odata, x1, x0, x2, x5, x4, x3);
}

void Task449::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I434
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, c2, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, c2, x0, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2, x0, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c2, x0, x2), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: l2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, x1)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
    dgemm_("T", "N", x2.size(), x1.size()*c2.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x2.size(), x0.size(), c2.size(), x1.size());
  out()->add_block(odata, x1, c2, x0, x2);
}

#endif
