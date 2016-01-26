//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks4.cc
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

#include <src/smith/relmrci/RelMRCI_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task150::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I267
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, c1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, c1, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c3, x5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c3, x5)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<2,1,0,3,0,1,1,2>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      zgemm3m_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, c1, x5);
}

void Task151::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          for (auto& x3 : *range_[1]) {
            // tensor label: Gamma92
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x5, x4, x2, x3, x1, x0);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x5, x4, x2, x3, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
            // tensor label: I285
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x3, c1, x7, x6);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x3, c1, x7, x6)]);
            sort_indices<4,5,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c1.size(), x7.size(), x6.size());
            zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*x7.size()*x6.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task152::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c1 = b(3);
  const Index x7 = b(4);
  const Index x6 = b(5);
  // tensor label: I285
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x5, x4, x3, c1, x7, x6)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x5, x4, x3, c1, x7, x6), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, x7, x6);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, x7, x6)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x7.size(), x6.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, a2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,1,2>(i1data, i1data_sorted, x5.size(), x4.size(), a2.size(), x3.size());
    zgemm3m_("T", "N", c1.size()*x7.size()*x6.size(), x5.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x7.size()*x6.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x7.size(), x6.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x5, x4, x3, c1, x7, x6);
}

void Task153::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        for (auto& x5 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            // tensor label: Gamma95
            std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x3, x6, x5, x2, x4, x1, x0);
            std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x3, x6, x5, x2, x4, x1, x0)]);
            sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x3.size(), x6.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
            // tensor label: I294
            std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x4, x3, x7, x6, x5);
            std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x4, x3, x7, x6, x5)]);
            sort_indices<3,2,4,5,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x7.size(), x6.size(), x5.size());
            zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, i0data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x7.size()*x6.size()*x5.size(),
                   1.0, odata_sorted, x2.size()*x1.size()*x0.size());
          }
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task154::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x7 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  // tensor label: I294
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x4, x3, x7, x6, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x4, x3, x7, x6, x5), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, a2, x6, x5);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, a2, x6, x5)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x4, a2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x4, a2, x3)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x4.size(), a2.size(), x3.size());
    zgemm3m_("T", "N", x7.size()*x6.size()*x5.size(), c1.size()*x4.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x7.size()*x6.size()*x5.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x4, x3, x7, x6, x5);
}

void Task155::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma98
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x3, x2, x4, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I303
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x4, x3, x5, c1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x4, x3, x5, c1)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), x5.size(), c1.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task156::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index x5 = b(2);
  const Index c1 = b(3);
  // tensor label: I303
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x3, x5, c1)]);
  std::fill_n(odata.get(), out()->get_size(x4, x3, x5, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x4, x3, x5, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, x5, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, a2, c1, a3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, a2, c1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), a3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a3, x4, a2, x3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a3, x4, a2, x3)]);
      sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, a3.size(), x4.size(), a2.size(), x3.size());
      zgemm3m_("T", "N", x5.size()*c1.size(), x4.size()*x3.size(), a3.size()*a2.size(),
             1.0, i0data_sorted, a3.size()*a2.size(), i1data_sorted, a3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size());
  out()->add_block(odata, x4, x3, x5, c1);
}

void Task157::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma414
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x7, x6, c1, x5);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task158::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x9 : *range_[1]) {
    for (auto& x8 : *range_[1]) {
      for (auto& x7 : *range_[1]) {
        // tensor label: Gamma415
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x9, x8, x2, x7, x1, x0);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x9, x8, x2, x7, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x9.size(), x8.size(), x2.size(), x7.size(), x1.size(), x0.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x9, x8, c1, x7);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x9, x8, c1, x7)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x9.size(), x8.size(), c1.size(), x7.size());
        zgemm3m_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x9.size()*x8.size()*x7.size(),
               1.0, i0data_sorted, x9.size()*x8.size()*x7.size(), i1data_sorted, x9.size()*x8.size()*x7.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task159::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0, c1, a2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0, c1, a2), 0.0);
  {
    // tensor label: I24
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c3, x0, a2);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c3.size(), x0.size(), a2.size());
  }
  out()->add_block(odata, c3, x0, c1, a2);
}

void Task160::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I25
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c3, x1, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c3, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x1.size(), x0.size());
    zgemm3m_("T", "N", a2.size(), c1.size()*c3.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task161::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I25
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x3, c3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,-2,1>(i1data, i1data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, c3, x1, x0);
}

void Task162::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I28
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
  zgemm3m_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task163::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I28
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, x0)]);
  std::fill_n(odata.get(), out()->get_size(c3, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c3, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->add_block(odata, c3, x0);
}

void Task164::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I31
  std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x0);
  std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
  zgemm3m_("T", "N", c3.size()*a2.size(), c1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task165::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I31
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->add_block(odata, c1, x0);
}

void Task166::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I34
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, c3, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), c1.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a2.size(), c3.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task167::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c4.size());
    zgemm3m_("T", "N", a2.size()*c3.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c3.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task168::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size());
    zgemm3m_("T", "N", c3.size()*a2.size()*x1.size(), c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size()*a2.size()*x1.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task169::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), c4.size());
    zgemm3m_("T", "N", a2.size()*c1.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, a2.size()*c1.size()*x1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task170::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a4.size(), a2.size());
    zgemm3m_("T", "N", c3.size()*c1.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c3.size()*c1.size()*x1.size());
  }
  sort_indices<1,3,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task171::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c4);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c3.size(), c4.size());
    zgemm3m_("T", "N", c1.size()*a2.size()*x1.size(), c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size()*a2.size()*x1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task172::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, a2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a2.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*x1.size(), a2.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*x1.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task173::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a2, c4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a2, c4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      zgemm3m_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task174::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a2, c5, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a2, c5, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, c3, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, c3, c4)]);
      sort_indices<3,1,0,2,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c5.size(), c3.size(), c4.size());
      zgemm3m_("T", "N", a2.size()*x1.size(), c1.size()*c3.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a2.size()*x1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a2.size(), x1.size(), c1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task175::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a4, c3, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a4, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c3.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task176::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a5, c3, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a5, c3, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c3.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      zgemm3m_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task177::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c5, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c5, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c5.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task178::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a5, c4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, a5, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a5.size(), c4.size());
      zgemm3m_("T", "N", c3.size()*x1.size(), c1.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task179::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task180::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c4, a5, c1, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c4, a5, c1, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a5.size(), c1.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<3,2,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      zgemm3m_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task181::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, c5, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, c5, a4, a2)]);
      sort_indices<2,1,0,3,0,1,-1,1>(i1data, i1data_sorted, c3.size(), c5.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), c5.size()*a4.size(),
             1.0, i0data_sorted, c5.size()*a4.size(), i1data_sorted, c5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task182::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c4, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c4, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), x1.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, a5, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, a5, c4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), a5.size(), c4.size());
      zgemm3m_("T", "N", c1.size()*x1.size(), c3.size()*a2.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task183::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a2, c4, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I505
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, c1, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, c1, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c1.size(), c4.size());
      zgemm3m_("T", "N", c3.size()*a2.size(), x1.size()*c1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x1.size(), c1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task184::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I505
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x1, c1, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c1, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x1, c1, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c1.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c1, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c1, c4);
}

void Task185::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a5, c4, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a5.size(), c4.size(), a2.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, x1, a5, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, x1, a5, c4)]);
      sort_indices<2,3,0,1,0,1,2,1>(i1data, i1data_sorted, c1.size(), x1.size(), a5.size(), c4.size());
      zgemm3m_("T", "N", c3.size()*a2.size(), c1.size()*x1.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x1.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task186::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c4, a2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c4.size(), a2.size());
      // tensor label: I514
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task187::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I514
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x1, c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c3, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c3, c4);
}

void Task188::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a2, c4, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a2, c4, a5)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a5.size());
      // tensor label: I517
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, c3, c4);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, c3, c4)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), c3.size(), c4.size());
      zgemm3m_("T", "N", c1.size()*a2.size(), x1.size()*c3.size(), a5.size()*c4.size(),
             1.0, i0data_sorted, a5.size()*c4.size(), i1data_sorted, a5.size()*c4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), c3.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task189::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index c4 = b(3);
  // tensor label: I517
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a5, x1, c3, c4)]);
  std::fill_n(odata.get(), out()->get_size(a5, x1, c3, c4), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a5, x1, c3, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), x1.size(), c3.size(), c4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c3, x1, a5, c4);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c3.size(), x1.size(), a5.size(), c4.size());
  }
  out()->add_block(odata, a5, x1, c3, c4);
}

void Task190::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a5 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a5, c3, a4);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a5, c3, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a4.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<0,2,1,3,0,1,-2,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task191::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I34
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a5);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a5.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a5, x1, a4, a2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a5, x1, a4, a2)]);
      sort_indices<2,0,1,3,0,1,2,1>(i1data, i1data_sorted, a5.size(), x1.size(), a4.size(), a2.size());
      zgemm3m_("T", "N", c1.size()*c3.size(), x1.size()*a2.size(), a5.size()*a4.size(),
             1.0, i0data_sorted, a5.size()*a4.size(), i1data_sorted, a5.size()*a4.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), a2.size());
  out()->add_block(odata, c1, a2, c3, x1);
}

void Task192::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I52
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x0.size(), x1.size());
    zgemm3m_("T", "N", c1.size(), c3.size()*a2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task193::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I52
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c3, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c3, a2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c3, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a2.size());
  out()->add_block(odata, c3, a2, x0, x1);
}

void Task194::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I55
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x0, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
    zgemm3m_("T", "N", c3.size(), c1.size()*a2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task195::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I55
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, a2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, a2, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, a2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma9
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c1.size(), a2.size());
  out()->add_block(odata, c1, a2, x0, x1);
}

void Task196::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I24
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I58
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    zgemm3m_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->add_block(odata, c1, c3, x0, a2);
}

void Task197::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I58
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma11
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, a4.size(), x1.size());
    zgemm3m_("T", "N", x0.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task198::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I58
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma160
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x3, x2, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x3, x2, x1)]);
        sort_indices<1,2,3,0,0,1,-1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

void Task199::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I58
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata.get(), out()->get_size(a4, x0), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma9
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: v2
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x3, x2, a4, x1);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x3, x2, a4, x1)]);
        sort_indices<0,1,3,2,0,1,-1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a4.size(), x1.size());
        zgemm3m_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->add_block(odata, a4, x0);
}

#endif
