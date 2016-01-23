//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks17.cc
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

#include <src/smith/mrci/MRCI_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task800::Task_local::compute() {
  const Index a1 = b(0);
  const Index c4 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1391
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c4, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c4, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c4, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c4.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c4.size(), a1.size());
  out()->add_block(odata, a1, c4, x0, x1);
}

void Task801::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, a4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, a4, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), a4.size(), a1.size());
      // tensor label: I1394
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c2, x0, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task802::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1394
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a4.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a4.size(), c2.size());
  out()->add_block(odata, a4, c2, x0, x1);
}

void Task803::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1394
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a4, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a4, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), a4.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a4.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a4.size());
  out()->add_block(odata, a4, c2, x0, x1);
}

void Task804::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, a4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, a4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), a4.size(), a3.size());
      // tensor label: I1397
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c2, x1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), a4.size()*x1.size(),
             1.0, i0data_sorted, a4.size()*x1.size(), i1data_sorted, a4.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task805::Task_local::compute() {
  const Index a4 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1397
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, c2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1398
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a4.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a4.size(), c2.size());
  out()->add_block(odata, a4, c2, x1, x0);
}

void Task806::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I1398
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a4, c2, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a4, c2, x2), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a4.size(), c2.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a4, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), a4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x3, a4, c2, x2);
}

void Task807::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
      // tensor label: I1400
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task808::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1400
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x0);
}

void Task809::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1400
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x0);
}

void Task810::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I1403
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x4, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x4.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task811::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1403
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x5, x4, x0);
}

void Task812::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1403
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<2,3,4,0,1,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, a3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, a3)]);
        sort_indices<0,1,2,3,0,1,-1,2>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), a3.size());
        dgemm_("T", "N", x5.size()*x4.size()*x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x5, x4, x0);
}

void Task813::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x2.size(), x1.size());
      // tensor label: I1430
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task814::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1430
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a3.size());
  out()->add_block(odata, a3, x0, x2, x1);
}

void Task815::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x2.size(), x1.size());
      // tensor label: I1433
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task816::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1433
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
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
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task817::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), x1.size(), a1.size());
      // tensor label: I1436
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task818::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1436
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x2, x1, x0);
}

void Task819::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2, x1, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size(), x1.size(), a3.size());
      // tensor label: I1439
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x1, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x1.size()*x2.size(),
             1.0, i0data_sorted, x1.size()*x2.size(), i1data_sorted, x1.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task820::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I1439
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma31
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), a1.size());
  out()->add_block(odata, a1, x0, x1, x2);
}

void Task821::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1, c2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1442
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1, x2, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x2.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), x1.size()*x2.size(),
             1.0, i0data_sorted, x1.size()*x2.size(), i1data_sorted, x1.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task822::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1442
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma471
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x4, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x1, x4, x3, x2, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,2>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x1, x2, x0);
}

void Task823::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a3, c2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x2.size(), a3.size(), c2.size(), x1.size());
      // tensor label: I1445
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task824::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1445
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task825::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I1448
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task826::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1448
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,-1,2>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a3.size());
  out()->add_block(odata, a3, x0, x2, x1);
}

void Task827::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x1.size(), c2.size(), a3.size());
      // tensor label: I1451
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task828::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1451
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
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
  out()->add_block(odata, a1, x0, x2, x1);
}

void Task829::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, a4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), a4.size(), a1.size());
    // tensor label: I1454
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task830::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I1454
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,2,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task831::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, a4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), a4.size(), a3.size());
    // tensor label: I1457
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task832::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I1457
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,-1,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task833::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c5, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c5, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c5.size(), a1.size());
      // tensor label: I1472
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c2, c4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c2, c4, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c5.size(), c2.size(), c4.size(), x0.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task834::Task_local::compute() {
  const Index c5 = b(0);
  const Index c2 = b(1);
  const Index c4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1472
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, c2, c4)]);
    sort_indices<0,1,2,3,0,1,4,1>(i1data, i1data_sorted, x1.size(), c5.size(), c2.size(), c4.size());
    dgemm_("T", "N", x0.size(), c5.size()*c2.size()*c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c5.size(), c2.size(), c4.size());
  out()->add_block(odata, c5, c2, c4, x0);
}

void Task835::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, c5, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a1, c5, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a1.size(), c5.size(), a3.size());
      // tensor label: I1475
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c2, c4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c2, c4, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c5.size(), c2.size(), c4.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task836::Task_local::compute() {
  const Index c5 = b(0);
  const Index c2 = b(1);
  const Index c4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1475
  std::unique_ptr<double[]> odata(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, c2, c4)]);
    sort_indices<0,1,2,3,0,1,-2,1>(i1data, i1data_sorted, x1.size(), c5.size(), c2.size(), c4.size());
    dgemm_("T", "N", x0.size(), c5.size()*c2.size()*c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c5.size(), c2.size(), c4.size());
  out()->add_block(odata, c5, c2, c4, x0);
}

void Task837::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), a1.size());
      // tensor label: I1502
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task838::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1502
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1503
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->add_block(odata, c2, c4, x3, x0);
}

void Task839::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1503
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, c2, c4);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x2.size(), x1.size(), c2.size(), c4.size());
  }
  out()->add_block(odata, c2, c4, x2, x1);
}

void Task840::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1502
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, x1, c4)]);
      sort_indices<1,2,0,3,0,1,1,2>(i1data, i1data_sorted, c2.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->add_block(odata, c2, c4, x3, x0);
}

void Task841::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1502
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c2, x1)]);
      sort_indices<3,0,1,2,0,1,-1,2>(i1data, i1data_sorted, x2.size(), c4.size(), c2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c4.size()*c2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c4.size(), c2.size());
  out()->add_block(odata, c2, c4, x3, x0);
}

void Task842::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c4.size(), a3.size());
      // tensor label: I1505
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task843::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1505
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1506
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->add_block(odata, c2, c4, x3, x0);
}

void Task844::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1506
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, c4, c2, x1);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x2.size(), c4.size(), c2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, c2, c4);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x2.size(), x1.size(), c2.size(), c4.size());
  }
  out()->add_block(odata, c2, c4, x2, x1);
}

void Task845::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1505
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,-1,2>(i1data, i1data_sorted, c2.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->add_block(odata, c2, c4, x3, x0);
}

void Task846::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
      // tensor label: I1508
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->add_block(odata, a3, c2, x0, a1);
}

void Task847::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1508
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1509
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

void Task848::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1509
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x2, x1), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), a1.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a1);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x2.size(), x1.size(), a4.size(), a1.size());
  }
  out()->add_block(odata, a4, a1, x2, x1);
}

void Task849::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1508
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->add_block(odata, a4, a1, x3, x0);
}

#endif
