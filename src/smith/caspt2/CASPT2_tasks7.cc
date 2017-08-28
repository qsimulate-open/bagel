//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_tasks7.cc
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

#include <src/smith/caspt2/CASPT2_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task300::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3), 0.0);
  {
    // tensor label: I372
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->add_block(odata, c2, x3);
}

void Task301::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I372
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
        // tensor label: I373
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
  out()->add_block(odata, x3, c2);
}

void Task302::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  // tensor label: I373
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x3, x4)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x3, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x3, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x3, x4), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x3, x4, x1, x0)]);
        sort_indices<0,4,5,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,2,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x3.size()*x4.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x3.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x3.size(), x4.size(), c1.size());
  out()->add_block(odata, c1, x5, x3, x4);
}

void Task303::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I372
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
        // tensor label: I529
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
  out()->add_block(odata, x3, c2);
}

void Task304::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  // tensor label: I529
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x3, x4)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x3, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x3, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x3, x4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma56
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x3.size()*x4.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x3.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x3.size(), x4.size(), a1.size());
  out()->add_block(odata, a1, x5, x3, x4);
}

void Task305::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I372
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
        // tensor label: I532
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
  out()->add_block(odata, x3, c2);
}

void Task306::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I532
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma57
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task307::Task_local::compute() {
  const Index x3 = b(0);
  const Index x4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, x4), 0.0);
  {
    // tensor label: I375
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x4.size(), x3.size());
  }
  out()->add_block(odata, x3, x4);
}

void Task308::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I375
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: Gamma143
              std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
              sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
              // tensor label: I376
              std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, x5, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, x5, x2, x1, x0)]);
              sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", x4.size()*x3.size(), 1, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(),
                     1.0, odata_sorted, x4.size()*x3.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->add_block(odata, x4, x3);
}

void Task309::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x5 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I376
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x7, x6, c1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x7.size(), x6.size(), x5.size());
  out()->add_block(odata, x7, x6, x5, x2, x1, x0);
}

void Task310::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  // tensor label: I375
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x2 : *range_[1]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: Gamma196
              std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
              sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
              // tensor label: I535
              std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, x5, x2, x1, x0);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, x5, x2, x1, x0)]);
              sort_indices<0,5,1,2,3,4,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x5.size(), x2.size(), x1.size(), x0.size());
              dgemm_("T", "N", x4.size()*x3.size(), 1, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(),
                     1.0, i0data_sorted, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(), i1data_sorted, x7.size()*x6.size()*x5.size()*x2.size()*x1.size()*x0.size(),
                     1.0, odata_sorted, x4.size()*x3.size());
            }
          }
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x4.size(), x3.size());
  out()->add_block(odata, x4, x3);
}

void Task311::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x5 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I535
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x7, x6, x5, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x7, x6, x5, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x7, a1, x6, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(0)->get_size(x7, a1, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x7.size(), x6.size(), x5.size());
  out()->add_block(odata, x7, x6, x5, x2, x1, x0);
}

void Task312::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1), 0.0);
  {
    // tensor label: I378
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c1.size(), c2.size());
  }
  out()->add_block(odata, c2, c1);
}

void Task313::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  // tensor label: I378
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I379
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
  out()->add_block(odata, c1, c2);
}

void Task314::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I379
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task315::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3), 0.0);
  {
    // tensor label: I381
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), a3.size());
  }
  out()->add_block(odata, c2, a3);
}

void Task316::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I381
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
      // tensor label: I382
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
      dgemm_("T", "N", c2.size()*a3.size(), 1, c1.size()*x3.size(),
             1.0, i0data_sorted, c1.size()*x3.size(), i1data_sorted, c1.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->add_block(odata, c2, a3);
}

void Task317::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I382
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma7
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,2,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->add_block(odata, c1, x3);
}

void Task318::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I381
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
      // tensor label: I385
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
      dgemm_("T", "N", a3.size()*c2.size(), 1, c1.size()*x3.size(),
             1.0, i0data_sorted, c1.size()*x3.size(), i1data_sorted, c1.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size());
  out()->add_block(odata, c2, a3);
}

void Task319::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I385
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma7
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size());
  out()->add_block(odata, c1, x3);
}

void Task320::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I381
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I541
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size());
      dgemm_("T", "N", a3.size()*c2.size(), 1, a1.size()*x3.size(),
             1.0, i0data_sorted, a1.size()*x3.size(), i1data_sorted, a1.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size());
  out()->add_block(odata, c2, a3);
}

void Task321::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  // tensor label: I541
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,-1,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, a1, x3);
}

void Task322::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  // tensor label: I381
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata.get(), out()->get_size(c2, a3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I544
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size());
      dgemm_("T", "N", c2.size()*a3.size(), 1, a1.size()*x3.size(),
             1.0, i0data_sorted, a1.size()*x3.size(), i1data_sorted, a1.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size());
  out()->add_block(odata, c2, a3);
}

void Task323::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  // tensor label: I544
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->add_block(odata, a1, x3);
}

void Task324::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  {
    // tensor label: I387
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x3.size(), a2.size());
  }
  out()->add_block(odata, x3, a2);
}

void Task325::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
        // tensor label: I388
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x3, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x3, x4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x3.size(), x4.size());
        dgemm_("T", "N", a2.size(), x3.size(), c1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c1.size()*x5.size()*x4.size(), i1data_sorted, c1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task326::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  // tensor label: I388
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x3, x4)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x3, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x3, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x3, x4), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x2, x4, x1, x0)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x3.size()*x4.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x3.size()*x4.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x3.size(), x4.size(), c1.size());
  out()->add_block(odata, c1, x5, x3, x4);
}

void Task327::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
        // tensor label: I391
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", a2.size(), x3.size(), c1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c1.size()*x5.size()*x4.size(), i1data_sorted, c1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task328::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma6
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<2,4,5,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task329::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  // tensor label: I387
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata.get(), out()->get_size(x3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
        // tensor label: I547
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", a2.size(), x3.size(), a1.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, a1.size()*x5.size()*x4.size(), i1data_sorted, a1.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size());
  out()->add_block(odata, x3, a2);
}

void Task330::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I547
  std::unique_ptr<double[]> odata(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma59
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<1,4,5,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1, x1, x2)]);
        sort_indices<0,3,2,1,0,1,2,1>(i1data, i1data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x5.size()*x4.size()*x3.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->add_block(odata, a1, x5, x4, x3);
}

void Task331::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata.get(), out()->get_size(a2, x1), 0.0);
  {
    // tensor label: I393
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), a2.size());
  }
  out()->add_block(odata, a2, x1);
}

void Task332::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  // tensor label: I393
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a2)]);
  std::fill_n(odata.get(), out()->get_size(x1, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I394
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x1, x0)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), x1.size(), c1.size()*c3.size()*x0.size(),
               1.0, i0data_sorted, c1.size()*c3.size()*x0.size(), i1data_sorted, c1.size()*c3.size()*x0.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size());
  out()->add_block(odata, x1, a2);
}

void Task333::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I394
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma3
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x1, x0);
}

void Task334::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  {
    // tensor label: I396
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), c1.size());
  }
  out()->add_block(odata, a2, c1);
}

void Task335::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I396
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata.get(), out()->get_size(a2, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I397
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), 1, c3.size()*x0.size(),
             1.0, i0data_sorted, c3.size()*x0.size(), i1data_sorted, c3.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->add_block(odata, a2, c1);
}

void Task336::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I397
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma12
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,2,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        dgemm_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task337::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(a2, c3)]);
  std::fill_n(odata.get(), out()->get_size(a2, c3), 0.0);
  {
    // tensor label: I399
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), a2.size());
  }
  out()->add_block(odata, a2, c3);
}

void Task338::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I399
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
      // tensor label: I400
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, c1.size()*x0.size(),
             1.0, i0data_sorted, c1.size()*x0.size(), i1data_sorted, c1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size());
  out()->add_block(odata, c3, a2);
}

void Task339::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I400
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma12
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task340::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I399
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I556
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size());
  out()->add_block(odata, c3, a2);
}

void Task341::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I556
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata.get(), out()->get_size(a4, c1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma38
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I557
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->add_block(odata, a4, c1);
}

void Task342::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I557
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, a4, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, a4, c1, x0), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x1.size(), a4.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x1, x0);
    sort_indices<2,1,0,3,1,1,-4,1>(i1data, odata, c1.size(), a4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x1, a4, c1, x0);
}

void Task343::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x2), 0.0);
  {
    // tensor label: I402
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
  }
  out()->add_block(odata, x1, x2);
}

void Task344::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I402
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma152
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
      // tensor label: I403
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
      dgemm_("T", "N", x2.size()*x1.size(), 1, x3.size()*x0.size(),
             1.0, i0data_sorted, x3.size()*x0.size(), i1data_sorted, x3.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->add_block(odata, x2, x1);
}

void Task345::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I403
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x0), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I404
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x3, x0);
}

void Task346::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I404
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, a2, c1, x3)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, c1, x3), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->add_block(odata, c3, a2, c1, x3);
}

void Task347::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I402
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I612
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
      dgemm_("T", "N", x2.size()*x1.size(), 1, x3.size()*x0.size(),
             1.0, i0data_sorted, x3.size()*x0.size(), i1data_sorted, x3.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->add_block(odata, x2, x1);
}

void Task348::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I612
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I613
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
        dgemm_("T", "N", x0.size(), x3.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->add_block(odata, x3, x0);
}

void Task349::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I613
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, a3, c2, a1)]);
  std::fill_n(odata.get(), out()->get_size(x3, a3, c2, a1), 0.0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->add_block(odata, x3, a3, c2, a1);
}

#endif
