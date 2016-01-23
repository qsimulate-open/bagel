//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks10.cc
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

#include <src/smith/caspt2/CASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task450::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I566
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I567
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task451::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I567
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
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
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

void Task452::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3), 0.0);
  {
    // tensor label: I572
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->add_block(odata, a4, c3);
}

void Task453::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I572
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I573
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
      dgemm_("T", "N", a4.size()*c3.size(), 1, a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, a4.size()*c3.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size());
  out()->add_block(odata, a4, c3);
}

void Task454::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I573
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma38
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I574
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->add_block(odata, a2, c1);
}

void Task455::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I574
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x1.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,1,0,3,1,1,8,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a2, c1, x0);
}

void Task456::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  {
    // tensor label: I581
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->add_block(odata, x0, x1);
}

void Task457::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I581
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  // tensor label: Gamma38
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I582
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", x1.size()*x0.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->add_block(odata, x1, x0);
}

void Task458::Task_local::compute() {
  // tensor label: I582
  std::unique_ptr<double[]> odata(new double[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: t2
  std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,3,2,1,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  odata_sorted[0] += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->add_block(odata);
}

void Task459::Task_local::compute() {
  // tensor label: I582
  std::unique_ptr<double[]> odata(new double[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: t2
  std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  odata_sorted[0] += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->add_block(odata);
}

void Task460::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c3)]);
  std::fill_n(odata.get(), out()->get_size(c5, c3), 0.0);
  {
    // tensor label: I587
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), c5.size());
  }
  out()->add_block(odata, c5, c3);
}

void Task461::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I587
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, a2)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,3,1,2,0,1,8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", c5.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, c5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c5.size(), c3.size());
  out()->add_block(odata, c3, c5);
}

void Task462::Task_local::compute() {
  const Index c3 = b(0);
  const Index c5 = b(1);
  // tensor label: I587
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata.get(), out()->get_size(c3, c5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, c5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, c5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,3,2,0,1,-16,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", c5.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, c5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c5.size(), c3.size());
  out()->add_block(odata, c3, c5);
}

void Task463::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata.get(), out()->get_size(a4, a5), 0.0);
  {
    // tensor label: I591
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a5);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a4.size(), a5.size());
  }
  out()->add_block(odata, a4, a5);
}

void Task464::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: I591
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata.get(), out()->get_size(a4, a5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,2,1,3,0,1,-8,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", a5.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a4, a5);
}

void Task465::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: I591
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata.get(), out()->get_size(a4, a5), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a5)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,16,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", a5.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->add_block(odata, a4, a5);
}

void Task466::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata.get(), out()->get_size(x0, c3), 0.0);
  {
    // tensor label: I595
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->add_block(odata, x0, c3);
}

void Task467::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I595
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I596
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task468::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I596
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<3,0,1,2,0,1,-4,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task469::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I596
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<1,0,3,2,0,1,2,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        dgemm_("T", "N", x1.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->add_block(odata, c3, x1);
}

void Task470::Task_local::compute() {
  const Index a1 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x1), 0.0);
  {
    // tensor label: I601
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a1.size());
  }
  out()->add_block(odata, a1, x1);
}

void Task471::Task_local::compute() {
  const Index x1 = b(0);
  const Index a1 = b(1);
  // tensor label: I601
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a1)]);
  std::fill_n(odata.get(), out()->get_size(x1, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I602
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x1, x0)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x1.size(), x0.size());
        dgemm_("T", "N", a1.size(), x1.size(), a3.size()*c2.size()*x0.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*x0.size(), i1data_sorted, a3.size()*c2.size()*x0.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x1.size());
  out()->add_block(odata, x1, a1);
}

void Task472::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I602
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I603
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), c2.size());
  out()->add_block(odata, a3, c2, x1, x0);
}

void Task473::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I603
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a3, c2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c2, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a3, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), a3.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a3, c2, x2);
}

void Task474::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1), 0.0);
  {
    // tensor label: I604
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a3.size());
  }
  out()->add_block(odata, a3, x1);
}

void Task475::Task_local::compute() {
  const Index x1 = b(0);
  const Index a3 = b(1);
  // tensor label: I604
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a3)]);
  std::fill_n(odata.get(), out()->get_size(x1, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I605
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c2, x0, x1)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), c2.size(), x0.size(), x1.size());
        dgemm_("T", "N", a3.size(), x1.size(), a1.size()*c2.size()*x0.size(),
               1.0, i0data_sorted, a1.size()*c2.size()*x0.size(), i1data_sorted, a1.size()*c2.size()*x0.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size());
  out()->add_block(odata, x1, a3);
}

void Task476::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I605
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task477::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I605
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->add_block(odata, a1, c2, x0, x1);
}

void Task478::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a1, c2), 0.0);
  {
    // tensor label: I613
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c2.size(), a1.size());
  }
  out()->add_block(odata, a1, c2);
}

void Task479::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  // tensor label: I613
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I614
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), 1, a3.size()*x0.size(),
             1.0, i0data_sorted, a3.size()*x0.size(), i1data_sorted, a3.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size());
  out()->add_block(odata, c2, a1);
}

void Task480::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  // tensor label: I614
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size());
  out()->add_block(odata, a3, x0);
}

void Task481::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2), 0.0);
  {
    // tensor label: I616
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a3.size(), c2.size());
  }
  out()->add_block(odata, a3, c2);
}

void Task482::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I616
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I617
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), 1, a1.size()*x0.size(),
             1.0, i0data_sorted, a1.size()*x0.size(), i1data_sorted, a1.size()*x0.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->add_block(odata, a3, c2);
}

void Task483::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  // tensor label: I617
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,2,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->add_block(odata, a1, x0);
}

void Task484::Task_local::compute() {
  const Index c4 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c4, x1), 0.0);
  {
    // tensor label: I619
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, x1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), x1.size());
  }
  out()->add_block(odata, c4, x1);
}

void Task485::Task_local::compute() {
  const Index c4 = b(0);
  const Index x1 = b(1);
  // tensor label: I619
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, x1)]);
  std::fill_n(odata.get(), out()->get_size(c4, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I620
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size());
    dgemm_("T", "N", x1.size(), c4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c4.size());
  out()->add_block(odata, c4, x1);
}

void Task486::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  // tensor label: I620
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<2,3,1,0,0,1,-4,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", c4.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size());
  out()->add_block(odata, x0, c4);
}

void Task487::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  // tensor label: I620
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata.get(), out()->get_size(x0, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<2,1,3,0,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        dgemm_("T", "N", c4.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size());
  out()->add_block(odata, x0, c4);
}

void Task488::Task_local::compute() {
  const Index c4 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c4, c2)]);
  std::fill_n(odata.get(), out()->get_size(c4, c2), 0.0);
  {
    // tensor label: I631
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->add_block(odata, c4, c2);
}

void Task489::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I631
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
        // tensor label: I632
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
  out()->add_block(odata, c2, c4);
}

void Task490::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I632
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->add_block(odata, a3, c2, a1, x1);
}

void Task491::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I631
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I635
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
  out()->add_block(odata, c2, c4);
}

void Task492::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I635
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,-2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->add_block(odata, a3, c2, a1, x1);
}

void Task493::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4), 0.0);
  {
    // tensor label: I637
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a1.size(), a4.size());
  }
  out()->add_block(odata, a1, a4);
}

void Task494::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: I637
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
        // tensor label: I638
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
  out()->add_block(odata, a1, a4);
}

void Task495::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I638
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->add_block(odata, a3, c2, a1, x1);
}

void Task496::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: I637
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata.get(), out()->get_size(a1, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
        // tensor label: I641
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
  out()->add_block(odata, a1, a4);
}

void Task497::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I641
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->add_block(odata, a3, c2, a1, x1);
}

void Task498::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a4)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4), 0.0);
  {
    // tensor label: I643
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a3.size(), a4.size());
  }
  out()->add_block(odata, a3, a4);
}

void Task499::Task_local::compute() {
  const Index a3 = b(0);
  const Index a4 = b(1);
  // tensor label: I643
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, a4)]);
  std::fill_n(odata.get(), out()->get_size(a3, a4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
        // tensor label: I644
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
  out()->add_block(odata, a3, a4);
}

#endif
