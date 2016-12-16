//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks3.cc
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

#include <src/smith/mrci/MRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task100::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma533
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x0, x6, x5, x4, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x5, x4, x1, x3, x2), 0.0);
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x2, x4, x1);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i5+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2+x2.size()*(i4+x4.size()*(i1)))))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i1+x1.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i3+x3.size()*(i2)))))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x7, x0, x6, x5, x4, x1, x3, x2);
}

void Task101::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma545
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x0, x4, x2, x3, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x2, x3, x1), 0.0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
  }
  out()->add_block(odata, x5, x0, x4, x2, x3, x1);
}

void Task102::Task_local::compute() {
  const Index x4 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma548
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x5, x1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5, x1, x4), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(4)->get_block(x3, x2);
  // rdm0 merged case
  if (x3 == x5 && x1 == x4 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          odata[i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4)))]  += -4.0 * i0data[0] * fdata[i5+x3.size()*(i2)];
        }
      }
    }
  }
  // rdm0 merged case
  if (x3 == x4 && x0 == x2 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          odata[i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4)))]  += 2.0 * i0data[0] * fdata[i4+x3.size()*(i2)];
        }
      }
    }
  }
  if (x0 == x2 && x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (-1.0) * i0data[i1+x1.size()*(i5)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x5 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (2.0) * i0data[i1+x1.size()*(i4)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x1 == x2 && x0 == x4 && x3 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          odata[i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]  += 2.0 * i0data[0] * fdata[i5+x3.size()*(i2)];
        }
      }
    }
  }
  if (x3 == x5 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (-1.0) * i0data[i1+x1.size()*(i2)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x1 == x2 && x3 == x4 && x0 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          odata[i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]  += -4.0 * i0data[0] * fdata[i4+x3.size()*(i2)];
        }
      }
    }
  }
  if (x3 == x4 && x0 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (2.0) * i0data[i1+x1.size()*(i2)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x1 == x2 && x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]
              += (2.0) * i0data[i0+x0.size()*(i5)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x5 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]
              += (-1.0) * i0data[i0+x0.size()*(i4)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x5 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4)))]
              += (2.0) * i0data[i0+x0.size()*(i2)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4)))]
              += (-1.0) * i0data[i0+x0.size()*(i2)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
                += (-1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i2)))] * fdata[i4+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  if (x3 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
                += (-1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i0+x0.size()*(i2)))] * fdata[i5+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x0, x5, x1, x4);
}

void Task103::Task_local::compute() {
  const Index x6 = b(0);
  const Index x1 = b(1);
  const Index x7 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  const Index x4 = b(6);
  const Index x5 = b(7);
  // tensor label: Gamma549
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x7, x1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x0, x7, x1, x6), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(5)->get_block(x5, x4, x3, x2);
  if (x0 == x2 && x3 == x7 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (-4.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x6 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x7 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x2 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x4 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i6+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x2 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-4.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x4 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i0+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x3 == x6 && x0 == x2 && x1 == x4 && x5 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]  += 2.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x4 && x0 == x2 && x3 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i7)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x1 == x4 && x0 == x2 && x5 == x6 && x3 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]  += -4.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x4 && x0 == x2 && x3 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i6)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x1 == x4 && x5 == x6 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x0 == x2 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x3 == x4 && x5 == x7 && x1 == x6 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]  += -4.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x6 && x5 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x6 && x3 == x4 && x1 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]  += 2.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x0 == x2 && x5 == x6 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (-1.0) * i0data[i1+x1.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (2.0) * i0data[i1+x1.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (-1.0) * i0data[i1+x1.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x7, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i7+x7.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x0 == x2 && x5 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (2.0) * i0data[i1+x1.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x6, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (2.0) * i0data[i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x7, x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i3+x3.size()*(i7+x7.size()*(i1+x1.size()*(i4)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x6, x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i3+x3.size()*(i4)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x7 && x3 == x6 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += -4.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i5+x5.size()*(i7)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x6 && x0 == x4 && x3 == x7 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += 2.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x0 == x4 && x3 == x7 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i6)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x6 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (-4.0) * i0data[i3+x3.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x7 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (2.0) * i0data[i1+x1.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x7, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i1+x1.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x7 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (-1.0) * i0data[i1+x1.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x6, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x7, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i7+x7.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x6, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (2.0) * i0data[i1+x1.size()*(i6+x6.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x7 && x3 == x4 && x1 == x2 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i6+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += 2.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x2 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x4 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (-1.0) * i0data[i1+x1.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x0 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i6+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x6 && x3 == x4 && x1 == x2 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += -4.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i3+x3.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x4 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-4.0) * i0data[i3+x3.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (2.0) * i0data[i1+x1.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x0 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (2.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i0+x0.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i0+x0.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (2.0) * i0data[i0+x0.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x7, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                  += (2.0) * i0data[i0+x0.size()*(i7+x7.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x7 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i0+x0.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x6, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x7, x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i3+x3.size()*(i4)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x6, x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                  += (-1.0) * i0data[i3+x3.size()*(i6+x6.size()*(i0+x0.size()*(i4)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i0+x0.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x7, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i5+x5.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x7 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (2.0) * i0data[i0+x0.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x6, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                  += (-1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i0+x0.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x7, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                  += (2.0) * i0data[i0+x0.size()*(i7+x7.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x6, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (2.0) * i0data[i0+x0.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                  += (2.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i0+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (-1.0) * i0data[i0+x0.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i0+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x7, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x6, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i0+x0.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x0, x7, x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i0+x0.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x7) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x1, x6, x5, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x0, x7, x1, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x1, x6, x0, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i0+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x0, x7, x1, x6);
}

void Task104::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // scalar
  // tensor label: Gamma558
  std::unique_ptr<double[]> odata(new double[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  out()->add_block(odata);
}

void Task105::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // scalar
  // tensor label: Gamma560
  std::unique_ptr<double[]> odata(new double[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(2)->get_block(x3, x2, x1, x0);
  out()->add_block(odata);
}

void Task106::Task_local::compute() {
  const Index x1 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma566
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x0, x4, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x1), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(2)->get_block(x3, x2);
  if (x4 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x5, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i5+x5.size()*(i0+x0.size()*(i2+x4.size()*(i1)))]
                += (1.0) * i0data[i3+x3.size()*(i1+x1.size()*(i5+x5.size()*(i0)))] * fdata[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  if (x5 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x4, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))]
                += (1.0) * i0data[i3+x3.size()*(i0+x0.size()*(i4+x4.size()*(i1)))] * fdata[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x5, x0, x4, x1);
}

void Task107::Task_local::compute() {
  const Index x1 = b(0);
  const Index x6 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  const Index x4 = b(6);
  const Index x5 = b(7);
  // tensor label: Gamma567
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x0, x6, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x1), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x5, x4, x3, x2);
  if (x7 == x4 && x6 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i4+x7.size()*(i0+x0.size()*(i2+x6.size()*(i1)))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i1)))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x2 && x6 == x4) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x3, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i2+x7.size()*(i0+x0.size()*(i4+x6.size()*(i1)))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i0)))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x3, x2, x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i0+x0.size()*(i4+x6.size()*(i1)))]
                    += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i2+x2.size()*(i7+x7.size()*(i0)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x3, x2, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i4+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))]
                    += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i1)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x6 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i0+x0.size()*(i2+x6.size()*(i1)))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i7+x7.size()*(i0)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x1, x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i7+x7.size()*(i0+x0.size()*(i2+x6.size()*(i1)))]
                    += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i7+x7.size()*(i0)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x7 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i6+x6.size()*(i1)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))]
                    += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i6+x6.size()*(i1)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x7, x0, x6, x1);
}

void Task109::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c1, x0, x1);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, c2.size(), c1.size(), x0.size(), x1.size());
  }
  {
    // tensor label: I0
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2, x1, x0);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task110::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I1
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task111::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c2.size(), c3.size());
    dgemm_("T", "N", c1.size()*x3.size()*x2.size(), c2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task112::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x3.size());
      // tensor label: I280
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task113::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I280
  std::unique_ptr<double[]> odata(new double[out()->get_size(a4, x2, c2, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c2, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), c2.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c2, c3);
}

void Task114::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x2.size(), a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), c2.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), c2.size(), x2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task115::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I4
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0, x1, x2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", c2.size(), c1.size()*x0.size()*x1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), c1.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task116::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I4
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, x0, x1, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0, x1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0, x1, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma1
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x0, x3, x1, x2)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        dgemm_("T", "N", x0.size()*x1.size()*x2.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x1.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), c1.size());
  out()->add_block(odata, c1, x0, x1, x2);
}

void Task117::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I7
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c1, c2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c1, c2, x3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c1.size(), c2.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*c2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task118::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I7
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, c1, c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, c1, c2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c1, c2, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,-1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", c1.size()*c2.size()*x3.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*c2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, x2, c1, c2, x3);
}

void Task119::Task_local::compute() {
  const Index x2 = b(0);
  const Index c1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I7
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, c1, c2, x3)]);
  std::fill_n(odata.get(), out()->get_size(x2, c1, c2, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c1, c2, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, x2, c1, c2, x3);
}

void Task120::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma80
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5, x1, x4, x3, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5, x1, x4, x3, x2)]);
          sort_indices<1,3,4,5,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
          // tensor label: I249
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c1, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c1, x5, x4)]);
          sort_indices<4,5,1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c1.size(), x5.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task121::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I249
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c3, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
    // tensor label: I250
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3, x3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x3.size(), x2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), c2.size()*x3.size()*x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, c1, x5, x4);
}

void Task122::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I250
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, c2, c3);
    sort_indices<2,3,0,1,1,1,-1,1>(i1data, odata, x3.size(), x2.size(), c2.size(), c3.size());
  }
  out()->add_block(odata, c2, c3, x3, x2);
}

void Task123::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma81
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5, x2, x4, x1, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5, x2, x4, x1, x3)]);
          sort_indices<1,2,3,5,0,4,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x2.size(), x4.size(), x1.size(), x3.size());
          // tensor label: I252
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c1, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c1, x5, x4)]);
          sort_indices<4,2,5,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c1.size(), x5.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task124::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I252
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c3, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c3)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c3.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), c2.size()*x3.size()*x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, c1, x5, x4);
}

void Task125::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma82
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5, x3, x4, x1, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5, x3, x4, x1, x2)]);
          sort_indices<1,2,3,5,0,4,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x3.size(), x4.size(), x1.size(), x2.size());
          // tensor label: I255
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2, x2, c1, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2, x2, c1, x5, x4)]);
          sort_indices<4,0,5,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size(), x2.size(), c1.size(), x5.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task126::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I255
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c2, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, c2, x2, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c2, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c2, x2, c1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c3, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, c2, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, c2, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), c2.size(), x2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size()*c2.size()*x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size(), c2.size(), x2.size());
  out()->add_block(odata, x3, c2, x2, c1, x5, x4);
}

void Task127::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        // tensor label: I264
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x6, x0, x5, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x6, x0, x5, x1)]);
        sort_indices<1,2,4,0,3,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x6.size(), x0.size(), x5.size(), x1.size());
        dgemm_("T", "N", c1.size(), c2.size()*x0.size()*x1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x0.size(), x1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task128::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I264
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x6, x0, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma85
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x0, x5, x1, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x0, x5, x1, x4, x3, x2)]);
        sort_indices<5,6,7,0,1,2,3,4,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x0.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,2>(i1data, i1data_sorted, c2.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x5.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x5, x1);
}

void Task129::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: I264
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x7, x6, x0, x5, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, x7, x6, x0, x5, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x5, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x5, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma86
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x0, x5, x4, x3, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x0, x5, x4, x3, x1, x2)]);
        sort_indices<4,5,7,0,1,2,3,6,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x0.size(), x5.size(), x4.size(), x3.size(), x1.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c2, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c2, x2)]);
        sort_indices<0,1,3,2,0,1,1,2>(i1data, i1data_sorted, x4.size(), x3.size(), c2.size(), x2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x5.size()*x1.size(), c2.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x5.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x5.size(), x1.size(), c2.size());
  out()->add_block(odata, c2, x7, x6, x0, x5, x1);
}

void Task130::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, c2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2, c2, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size(), c2.size(), c3.size());
      // tensor label: I270
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1, x0, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size(), x0.size(), x2.size());
      dgemm_("T", "N", c1.size()*c2.size(), x1.size()*x0.size(), c3.size()*x2.size(),
             1.0, i0data_sorted, c3.size()*x2.size(), i1data_sorted, c3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task131::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I270
  std::unique_ptr<double[]> odata(new double[out()->get_size(c3, x1, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(c3, x1, x0, x2), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1, x0, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1, x0, x2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma87
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x0, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x0, x2)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x0.size(), x2.size());
        // tensor label: t2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c3, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c3.size(), x3.size());
        dgemm_("T", "N", x1.size()*x0.size()*x2.size(), c3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x0.size()*x2.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), x2.size(), c3.size());
  out()->add_block(odata, c3, x1, x0, x2);
}

void Task132::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x5)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x5.size());
      // tensor label: I273
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x1, x5, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x1, x5, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x1.size(), x5.size(), x0.size());
      dgemm_("T", "N", c1.size()*c2.size(), x1.size()*x0.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x1.size(), x0.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task133::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  // tensor label: I273
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x1, x5, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1, x5, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x5, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x5, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma88
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x0, x4, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x0, x4, x3, x2)]);
        sort_indices<3,4,5,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,-1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x1.size()*x5.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x1.size()*x5.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x5.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x1, x5, x0);
}

void Task134::Task_local::compute() {
  const Index a3 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  // tensor label: I273
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x1, x5, x0)]);
  std::fill_n(odata.get(), out()->get_size(a3, x1, x5, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x1, x5, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x1, x5, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma89
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x5, x4, x3, x0, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x5, x4, x3, x0, x2)]);
        sort_indices<2,3,5,0,1,4,0,1,1,1>(i0data, i0data_sorted, x1.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size());
        // tensor label: v2
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<0,1,3,2,0,1,-1,2>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x1.size()*x5.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x1.size()*x5.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x5.size(), x0.size(), a3.size());
  out()->add_block(odata, a3, x1, x5, x0);
}

void Task135::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma94
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x0, x4, x1, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x0, x4, x1, x3)]);
          sort_indices<0,1,3,5,2,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x0.size(), x4.size(), x1.size(), x3.size());
          // tensor label: I291
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, x5, c1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, x5, c1, x4)]);
          sort_indices<3,2,5,1,0,4,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), x5.size(), c1.size(), x4.size());
          dgemm_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task136::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c1 = b(4);
  const Index x4 = b(5);
  // tensor label: I291
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x2, x5, c1, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, x5, c1, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c1.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), c2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, x5, c1, x4);
}

void Task137::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma87
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x0, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x0, x2)]);
          sort_indices<0,1,3,5,2,4,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x0.size(), x2.size());
          // tensor label: I294
          std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x2, c1, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x2, c1, x5, x4)]);
          sort_indices<4,5,1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c1.size(), x5.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task138::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I294
  std::unique_ptr<double[]> odata(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x5.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, a3, x2)]);
    sort_indices<2,0,1,3,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x3.size(), a3.size(), x2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), c2.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, c1, x5, x4);
}

void Task139::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x0, x1), 0.0);
  {
    // tensor label: I9
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c1, x2, x0, x1);
}

void Task140::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma3
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x5, x3, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
        // tensor label: I10
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c1, x5, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c1, x5, x4)]);
        sort_indices<2,0,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), c1.size(), x5.size(), x4.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task141::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I10
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,2,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, c1, x5, x4);
}

void Task142::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I10
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x5, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x5, c3, x4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x5.size(), c3.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, c1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, c1, c2)]);
      sort_indices<3,1,0,2,0,1,-2,1>(i1data, i1data_sorted, x3.size(), c3.size(), c1.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*c1.size(), c3.size()*c2.size(),
             1.0, i0data_sorted, c3.size()*c2.size(), i1data_sorted, c3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, x3, c1, x5, x4);
}

void Task143::Task_local::compute() {
  const Index x3 = b(0);
  const Index c1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I10
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, c1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x5)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x5.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<3,0,1,2,0,1,1,2>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", c1.size()*x5.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, c1.size()*x5.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, x3, c1, x5, x4);
}

void Task144::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I9
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma4
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I13
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->add_block(odata, c1, x2, x1, x0);
}

void Task145::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I13
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c2)]);
    sort_indices<1,0,0,1,-1,1>(i1data, i1data_sorted, c1.size(), c2.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), c1.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task146::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I13
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
    // tensor label: h1
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task147::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I13
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
      // tensor label: v2
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, a3, c2)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c1.size(), x3.size(), a3.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), c1.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), c1.size(), x3.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task148::Task_local::compute() {
  const Index c1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I13
  std::unique_ptr<double[]> odata(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
      // tensor label: I370
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, c1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, c1, c2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), c1.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*c1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), c1.size());
  out()->add_block(odata, c1, x5, x4, x3);
}

void Task149::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I370
  std::unique_ptr<double[]> odata(new double[out()->get_size(a3, x3, c1, c2)]);
  std::fill_n(odata.get(), out()->get_size(a3, x3, c1, c2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, c1, c2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a3.size(), x3.size(), c1.size(), c2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x3, a3, c2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c1.size(), x3.size(), a3.size(), c2.size());
  }
  out()->add_block(odata, a3, x3, c1, c2);
}

#endif
