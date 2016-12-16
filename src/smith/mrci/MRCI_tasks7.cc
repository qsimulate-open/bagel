//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks7.cc
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

#include <src/smith/mrci/MRCI_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task300::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), x4.size());
      // tensor label: I537
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x0, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x0, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x0.size(), x4.size());
      dgemm_("T", "N", a2.size()*c3.size(), c1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task301::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I537
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma176
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x0, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x0, x4, x2, x1)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x0.size(), x4.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c1.size());
  out()->add_block(odata, c1, x5, x0, x4);
}

void Task302::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  // tensor label: I537
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x0, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x0, x4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma178
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x0, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x0, x4, x3, x2)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x0.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x5.size()*x0.size()*x4.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x0.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x0.size(), x4.size(), c1.size());
  out()->add_block(odata, c1, x5, x0, x4);
}

void Task303::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I540
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task304::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I540
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,2>(i1data, i1data_sorted, c3.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x4, x0);
}

void Task305::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I540
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x4, x0);
}

void Task306::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), c4.size());
      // tensor label: I549
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c4, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c4, x0, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c4.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*c1.size(), a2.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task307::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I549
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I550
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c4, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c4.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c4.size());
  out()->add_block(odata, a2, c4, x0, x1);
}

void Task308::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c4 = b(2);
  const Index x2 = b(3);
  // tensor label: I550
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2, c4, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2, c4, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c4, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c4.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c4, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c4.size(), a2.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a2, c4, x2);
}

void Task309::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c1.size(), a2.size());
      // tensor label: I552
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x1, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task310::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I552
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c3, x2)]);
      sort_indices<0,3,1,2,0,1,-2,1>(i1data, i1data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c3.size());
  out()->add_block(odata, a4, c3, x1, x0);
}

void Task311::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I552
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,2,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a4.size());
  out()->add_block(odata, a4, c3, x1, x0);
}

void Task312::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), c4.size());
      // tensor label: I555
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c4, x1, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), c4.size()*x1.size(),
             1.0, i0data_sorted, c4.size()*x1.size(), i1data_sorted, c4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task313::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I555
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c4, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c4, x2)]);
      sort_indices<0,3,1,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c4.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c4.size());
  out()->add_block(odata, a2, c4, x1, x0);
}

void Task314::Task_local::compute() {
  const Index a2 = b(0);
  const Index c4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I555
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c4.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c4.size(), a2.size());
  out()->add_block(odata, a2, c4, x1, x0);
}

void Task315::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, a4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), a4.size(), a2.size());
      // tensor label: I558
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task316::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c3.size());
  out()->add_block(odata, a4, c3, x1, x0);
}

void Task317::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I558
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a4.size());
  out()->add_block(odata, a4, c3, x1, x0);
}

void Task318::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size(), c3.size(), a2.size());
      // tensor label: I561
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1, x1, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task319::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I561
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c1.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c1.size());
  out()->add_block(odata, a4, c1, x1, x0);
}

void Task320::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I561
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c1.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), a4.size());
  out()->add_block(odata, a4, c1, x1, x0);
}

void Task321::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, a4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), a4.size(), a2.size());
      // tensor label: I564
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1, x0, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task322::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I564
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I565
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a4.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a4.size(), c1.size());
  out()->add_block(odata, a4, c1, x0, x1);
}

void Task323::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I565
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a4, c1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a4, c1, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a4.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a4, c1, x2);
}

void Task324::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I567
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task325::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I567
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,2>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x0);
}

void Task326::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I567
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x0);
}

void Task327::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I570
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task328::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I570
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma132
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x4, x0);
}

void Task329::Task_local::compute() {
  const Index c3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I570
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma179
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x0, x1)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), c3.size());
  out()->add_block(odata, c3, x5, x4, x0);
}

void Task330::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2, c3, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size(), c3.size(), x1.size());
      // tensor label: I597
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*c3.size(), a2.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task331::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I597
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma196
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x0, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x0.size()*x1.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x0.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x0.size(), x1.size(), a2.size());
  out()->add_block(odata, a2, x2, x0, x1);
}

void Task332::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I642
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task333::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I642
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c1.size(), a4.size());
  out()->add_block(odata, c1, a4, x3, x0);
}

void Task334::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I645
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c3.size(), c1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task335::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I645
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, a4, x1)]);
      sort_indices<1,3,0,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c1.size(), a4.size());
  out()->add_block(odata, c1, a4, x3, x0);
}

void Task336::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a2.size());
      // tensor label: I648
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task337::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I648
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, c3.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x0);
}

void Task338::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a4.size());
      // tensor label: I651
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task339::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I651
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), a4.size());
  out()->add_block(odata, c3, a4, x3, x0);
}

void Task340::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma552
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size());
    // tensor label: I1681
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task341::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1681
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task342::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma554
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size());
    // tensor label: I1685
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task343::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1685
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x5), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x5);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), a2.size(), c1.size(), x5.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x5);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x5.size());
  }
  out()->add_block(odata, c3, a2, c1, x5);
}

void Task344::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  {
    // tensor label: I72
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, x0, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), x1.size(), x0.size(), a1.size());
  }
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task345::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I73
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, x2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, x2, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
    dgemm_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task346::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I73
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), c2.size());
  out()->add_block(odata, c2, x1, x2, x0);
}

void Task347::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I76
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, c2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, c2, x3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), c2.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->add_block(odata, c2, x1, x0, a1);
}

void Task348::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", a1.size()*c2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, a1, c2, x3);
}

void Task349::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a1, c3, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a1.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c2, c3)]);
      sort_indices<1,3,0,2,0,1,-1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c2.size(), c3.size());
      dgemm_("T", "N", a1.size()*x3.size(), x2.size()*c2.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a1.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, x2, a1, c2, x3);
}

#endif
