//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_contract_tasks.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_contract_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

void Task901::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  
  std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i0data, i0data_sorted, ci0.size());

  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", ci0.size(), 1, 1, 1.0, i0data_sorted, 1,
         i1data_sorted, 1, 1.0, odata_sorted, ci0.size());

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task902::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size());

      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(),
             i1data_sorted, x0.size()*x1.size(), 1.0, odata_sorted, ci0.size());
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task903::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size());

          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size(),
                 i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size(), 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}


void Task904::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3, x4, x5)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());

              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3, x4, x5)]);
              sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
              dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(), 1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}

void Task905::Task_local::compute() {
  const Index ci0 = b(0);

  std::unique_ptr<double[]> odata(new double[out()->get_size(ci0)]);
  std::fill_n(odata.get(), out()->get_size(ci0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1, x2, x3, x4, x5)]);
              sort_indices<1,2,3,4,5,6,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());

              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, x3, x4, x5);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, x3, x4, x5)]);
              sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), x3.size(), x4.size(), x5.size());
              dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     1.0, i0data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(),
                     i1data_sorted, x0.size()*x1.size()*x2.size()*x3.size()*x4.size()*x5.size(), 1.0, odata_sorted, ci0.size());
            }
          }
        }
      }
    }
  }

  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->add_block(odata, ci0);
}
#endif
