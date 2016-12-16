//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_tasks1.cc
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

#include <src/smith/mrci/MRCI_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task0::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: Gamma0
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x1, x2), 0.0);
  {
    // rdm0 non-merged case
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i2+x0.size()*(i3+x3.size()*(i3+x1.size()*(i2)))]  += -2.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i2+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))]
              += (1.0) * i0data[i1+x1.size()*(i3)];
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i3+x0.size()*(i3+x3.size()*(i2+x1.size()*(i2)))]  += 4.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))]
              += (-2.0) * i0data[i1+x1.size()*(i2)];
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            odata[i0+x0.size()*(i3+x3.size()*(i2+x1.size()*(i2)))]
              += (-2.0) * i0data[i0+x0.size()*(i3)];
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            odata[i0+x0.size()*(i3+x3.size()*(i3+x1.size()*(i2)))]
              += (1.0) * i0data[i0+x0.size()*(i2)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x3, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x0, x3, x1, x2);
}

void Task1::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma1
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x0, x3, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x0, x3, x1, x2), 0.0);
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i2+x0.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))]
                += (-2.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i2+x0.size()*(i3+x3.size()*(i4+x1.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x0.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))]
                += (4.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x0.size()*(i3+x3.size()*(i4+x1.size()*(i2)))))]
                += (-2.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x0.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))]
                += (-2.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x0.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i4+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x0, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x5, x4, x0, x3, x1, x2);
}

void Task2::Task_local::compute() {
  const Index x2 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  // tensor label: Gamma2
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x3, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x3, x0, x2), 0.0);
  {
    // rdm0 non-merged case
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i3+x1.size()*(i3+x3.size()*(i2+x0.size()*(i2)))]  += -4.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            odata[i1+x1.size()*(i3+x3.size()*(i2+x0.size()*(i2)))]
              += (2.0) * i0data[i1+x1.size()*(i3)];
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i2+x1.size()*(i3+x3.size()*(i3+x0.size()*(i2)))]  += 2.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            odata[i1+x1.size()*(i3+x3.size()*(i3+x0.size()*(i2)))]
              += (-1.0) * i0data[i1+x1.size()*(i2)];
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i2+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))]
              += (-1.0) * i0data[i0+x0.size()*(i3)];
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))]
              += (2.0) * i0data[i0+x0.size()*(i2)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x3, x0, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->add_block(odata, x1, x3, x0, x2);
}

void Task3::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x5 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma80
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x5, x1, x4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5, x1, x4, x3, x2), 0.0);
  {
    if (x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-2.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (4.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x4 && x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]  += 4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-2.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]  += -2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (1.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (1.0) * i0data[i1+x1.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (-2.0) * i0data[i1+x1.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x5 && x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]  += -2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (1.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (1.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]  += 4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-2.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (-2.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-2.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (-2.0) * i0data[i0+x0.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (1.0) * i0data[i0+x0.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (-2.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-2.0) * i0data[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (1.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                  += (1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                  += (1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x0, x5, x1, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x0.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x0, x5, x1, x4, x3, x2);
}

void Task4::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma81
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x5, x2, x4, x1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5, x2, x4, x1, x3), 0.0);
  {
    // rdm0 non-merged case
    if (x2 == x5 && x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i3+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]  += -2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (1.0) * i0data[i2+x2.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x5 && x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i3+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]  += 4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]
                += (-2.0) * i0data[i2+x2.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (-2.0) * i0data[i1+x1.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (1.0) * i0data[i1+x1.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x2, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x5 && x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]  += 4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (-2.0) * i0data[i2+x2.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x3 && x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]  += -2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]
                += (1.0) * i0data[i2+x2.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (1.0) * i0data[i1+x1.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (-2.0) * i0data[i1+x1.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x5, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x4 && x1 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]  += -8.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (4.0) * i0data[i2+x2.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x3 && x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]  += 4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (-2.0) * i0data[i2+x2.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (-2.0) * i0data[i1+x1.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                += (4.0) * i0data[i1+x1.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x4, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (-2.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (4.0) * i0data[i0+x0.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (-2.0) * i0data[i0+x0.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x2, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                  += (-2.0) * i0data[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (-2.0) * i0data[i0+x0.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (1.0) * i0data[i0+x0.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x2, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                  += (1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]
                += (1.0) * i0data[i0+x0.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]
                += (-2.0) * i0data[i0+x0.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x4, x0, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i3)))))]
                  += (1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x1, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (-2.0) * i0data[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x0, x5, x2, x4, x1, x3);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x0.size(), x5.size(), x2.size(), x4.size(), x1.size(), x3.size());
  }
  out()->add_block(odata, x0, x5, x2, x4, x1, x3);
}

void Task5::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma82
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x5, x3, x4, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5, x3, x4, x1, x2), 0.0);
  {
    if (x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i5+x1.size()*(i2)))))]
                += (-1.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x4 && x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]  += 2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]
                += (-1.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i5+x1.size()*(i2)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i5+x1.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                += (2.0) * i0data[i1+x1.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                += (-1.0) * i0data[i1+x1.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i2+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x5 && x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                += (2.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]  += 8.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                += (-4.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                += (-4.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x4, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (2.0) * i0data[i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                += (-4.0) * i0data[i0+x0.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                += (2.0) * i0data[i0+x0.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                  += (2.0) * i0data[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]
                += (-1.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              odata[i0+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i5+x1.size()*(i2)))))]
                += (2.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i5+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (2.0) * i0data[i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                odata[i0+x0.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x0, x5, x3, x4, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x0.size(), x5.size(), x3.size(), x4.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x0, x5, x3, x4, x1, x2);
}

void Task6::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x5 = b(4);
  const Index x0 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma85
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x0, x5, x1, x4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x0, x5, x1, x4, x3, x2), 0.0);
  {
    if (x1 == x4 && x0 == x2 && x3 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2 && x3 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i4)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i4)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x3, x5, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i4)))))];
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
    if (x3 == x5 && x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x6 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x1, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i4)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (4.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x6 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x1, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4 && x1 == x2 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i4)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x1, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i4)))))];
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
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i0+x0.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x4, x0, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i2)))))];
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
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x1, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i0+x0.size()*(i2)))))];
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
    if (x3 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x2, x0, x5, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i6+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x0, x5, x1, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x0.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x7, x6, x0, x5, x1, x4, x3, x2);
}

void Task7::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x0 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma86
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x0, x5, x4, x3, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x0, x5, x4, x3, x1, x2), 0.0);
  {
    if (x1 == x6 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2 && x1 == x3 && x4 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x2 && x4 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x1, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x5 && x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i6)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i5)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x6 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x4, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i2)))))];
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
    if (x4 == x6 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (4.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))];
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
    if (x4 == x5 && x1 == x2 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i2)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x0 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x4, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))];
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
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x5, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i2)))))];
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
    if (x4 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i0+x0.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))];
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
    if (x4 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x2, x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x2, x0, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x0, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))];
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
    if (x4 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x0, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i0+x0.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i2)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x0, x5, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x0.size(), x5.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x7, x6, x0, x5, x4, x3, x1, x2);
}

void Task8::Task_local::compute() {
  const Index x2 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma87
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x1, x3, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x1, x3, x0, x2), 0.0);
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                += (-4.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                += (2.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i2+x1.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                += (2.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                += (-1.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i2+x1.size()*(i3+x3.size()*(i4+x0.size()*(i2)))))]
                += (-1.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i4+x0.size()*(i2)))))]
                += (2.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i4+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x1, x3, x0, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->add_block(odata, x5, x4, x1, x3, x0, x2);
}

void Task9::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma88
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x5, x0, x4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x5, x0, x4, x3, x2), 0.0);
  {
    if (x1 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-4.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x1.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x1 == x4 && x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]  += 2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-1.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (2.0) * i0data[i1+x1.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (-1.0) * i0data[i1+x1.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i2+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x5 && x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x1.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (2.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (2.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (2.0) * i0data[i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x1.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]  += 2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                += (-1.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (-1.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (-1.0) * i0data[i0+x0.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (2.0) * i0data[i0+x0.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x0, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i2+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                += (-1.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (2.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (2.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x1, x5, x0, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x1.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x1, x5, x0, x4, x3, x2);
}

void Task10::Task_local::compute() {
  const Index x2 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma89
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x5, x4, x3, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x1, x5, x4, x3, x0, x2), 0.0);
  {
    if (x1 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                += (-4.0) * i0data[i4+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                += (2.0) * i0data[i4+x4.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x0.size()*(i2)))))]
                += (2.0) * i0data[i4+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x0.size()*(i2)))))]
                += (-1.0) * i0data[i4+x4.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x3, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (2.0) * i0data[i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x0 == x2 && x1 == x3 && x4 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i3+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                += (2.0) * i0data[i4+x4.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                += (2.0) * i0data[i1+x1.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x0.size()*(i2)))))]
                  += (2.0) * i0data[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x4 == x5 && x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]  += 2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                += (-1.0) * i0data[i4+x4.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                += (-1.0) * i0data[i1+x1.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x5, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                += (-1.0) * i0data[i0+x0.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                += (2.0) * i0data[i0+x0.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x5, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x1, x3, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i3+x3.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x1, x5, x4, x3, x0, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x1.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size());
  }
  out()->add_block(odata, x1, x5, x4, x3, x0, x2);
}

void Task11::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma94
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x2, x0, x4, x1, x3)]);
  std::fill_n(odata.get(), out()->get_size(x5, x2, x0, x4, x1, x3), 0.0);
  {
    if (x1 == x3 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i2+x0.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (-2.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i2+x0.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (1.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i2+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i3+x0.size()*(i4+x4.size()*(i2+x1.size()*(i3)))))]
                += (1.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i3+x0.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (-2.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x1, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i3+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i4+x0.size()*(i4+x4.size()*(i2+x1.size()*(i3)))))]
                += (-2.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i4+x0.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (4.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x0, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i4+x4.size()*(i2+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x0, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x0, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x2, x0, x4, x1, x3);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x2.size(), x0.size(), x4.size(), x1.size(), x3.size());
  }
  out()->add_block(odata, x5, x2, x0, x4, x1, x3);
}

void Task12::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma3
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x5, x3, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x5, x3, x4, x1, x0), 0.0);
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                += (-2.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                += (1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x2.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                += (1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                  += (1.0) * i0data[i2+x2.size()*(i5+x5.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x2.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                += (-2.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                  += (1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x2.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                += (-2.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x2.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                += (4.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x2.size()*(i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i5+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x5, x3, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x2, x5, x3, x4, x1, x0);
}

void Task13::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma4
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x2, x3, x1, x0), 0.0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                += (-1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))]
                += (2.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task14::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma5
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x1, x0), 0.0);
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))]
              += (-1.0) * i0data[i2+x2.size()*(i0)];
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))]
              += (2.0) * i0data[i1+x1.size()*(i0)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x2, x3, x1, x0);
}

void Task15::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma7
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x3, x2, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x3, x2, x4, x1, x0), 0.0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i4+x2.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))]
                += (-2.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i3+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x3, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x5, x3, x2, x4, x1, x0);
}

void Task16::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  const Index x7 = b(6);
  const Index x2 = b(7);
  // tensor label: Gamma97
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x7, x5, x6, x4, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x7, x5, x6, x4, x3, x1, x0), 0.0);
  {
    if (x2 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x4 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7 && x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x3 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x3 && x5 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i4+x4.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x2 == x7 && x5 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i4+x4.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x4 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i2+x2.size()*(i7+x7.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i0)))))];
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
    if (x5 == x7 && x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i4+x4.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x7, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i4+x4.size()*(i7+x7.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))))];
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
    if (x5 == x6 && x2 == x3 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i4+x4.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x6, x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i7+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))];
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
    if (x5 == x7 && x4 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x4 == x7 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i3+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i4+x4.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i4+x4.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x7, x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i3+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i4+x4.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x7, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i6+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (4.0) * i0data[i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7 && x4 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x4 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x6, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i7+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i2+x2.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i7+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x2, x7, x5, x6, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x2.size(), x7.size(), x5.size(), x6.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x2, x7, x5, x6, x4, x3, x1, x0);
}

void Task17::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x5 = b(3);
  const Index x6 = b(4);
  const Index x3 = b(5);
  const Index x7 = b(6);
  const Index x2 = b(7);
  // tensor label: Gamma98
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x7, x3, x6, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x7, x3, x6, x5, x4, x1, x0), 0.0);
  {
    if (x1 == x4 && x2 == x6 && x3 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x7 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7 && x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x4 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i5+x5.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (4.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x2 == x6 && x5 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i3+x3.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x7, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i7+x7.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4 && x2 == x7 && x5 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i3+x3.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x7 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x3 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x3 == x7 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x3, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i0)))))];
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
    if (x5 == x7 && x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i3+x3.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x3 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x4 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i3+x3.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x3 == x4 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i2+x2.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x6, x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))))];
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
    if (x5 == x7 && x3 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x3 == x7 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i4+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i3+x3.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x7, x3, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i4+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7 && x3 == x4 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i6+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i5+x5.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i6+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x7, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i6+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i3+x3.size()*(i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x3 == x4 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i1+x1.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x2 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x6, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x7, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x7 && x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7 && x3 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i2+x2.size()*(i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x3 == x7) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i2+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x6, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i7+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x7, x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x7) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x6, x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i2 = 0; i2 != x2.size(); ++i2) {
                    odata[i2+x2.size()*(i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i7+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i3+x3.size()*(i6+x6.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x2, x7, x3, x6, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x2.size(), x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x2, x7, x3, x6, x5, x4, x1, x0);
}

void Task18::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma100
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x2, x5, x4, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x2, x5, x4, x3, x1, x0), 0.0);
  {
    if (x1 == x5 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))];
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
    if (x4 == x5 && x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x2, x5, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x5 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x2, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x6, x2, x5, x4, x3, x1, x0);
}

void Task19::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x3 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma101
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x3, x5, x2, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x3, x5, x2, x4, x1, x0), 0.0);
  {
    if (x1 == x4 && x2 == x5 && x3 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x2 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x3, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x3, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))];
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
    if (x3 == x6 && x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))];
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
    if (x3 == x5 && x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (-4.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x3, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x3, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x3, x5, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4)))))];
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
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-4.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x3, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x3 == x4 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x4, x3, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i6+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x3 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x3, x5, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x6, x3, x5, x2, x4, x1, x0);
}

void Task20::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma102
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x4, x5, x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x4, x5, x2, x3, x1, x0), 0.0);
  {
    if (x2 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i6+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i3)))))];
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
    if (x4 == x5 && x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                  += (4.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x4, x5, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3)))))];
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
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (4.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x5 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i6+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x4, x5, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x6, x4, x5, x2, x3, x1, x0);
}

void Task21::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma104
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x5, x4, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x2, x5, x4, x3, x1, x0), 0.0);
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                += (2.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))]
                += (-1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (2.0) * i0data[i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                += (-1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                += (2.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x2.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x2, x5, x4, x3, x1, x0);
}

void Task22::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma107
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, x2, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, x2, x4, x1, x0), 0.0);
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                += (-1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                += (2.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                += (2.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                += (-1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                += (-4.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (2.0) * i0data[i3+x3.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                += (2.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (2.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x3.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x3, x5, x2, x4, x1, x0);
}

void Task23::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: Gamma109
  std::unique_ptr<double[]> odata(new double[out()->get_size(x4, x5, x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x4, x5, x2, x3, x1, x0), 0.0);
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                += (1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))]
                += (-2.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))]
                  += (1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i5+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                += (-2.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                += (4.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x4.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x5, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x4.size(), x5.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x4, x5, x2, x3, x1, x0);
}

void Task24::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x6 = b(4);
  const Index x2 = b(5);
  const Index x5 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma114
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x5, x2, x6, x4, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x5, x2, x6, x4, x3, x1, x0), 0.0);
  {
    if (x2 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i6+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i5+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3 && x4 == x6 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i5+x5.size()*(i5+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i5+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i5+x5.size()*(i6+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i6+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x6);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i6)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x2, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x4, x6);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i6)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x2, x6, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x5 && x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))];
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
    if (x4 == x5 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x4, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i3+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i5+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i5+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x5 && x2 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i5+x5.size()*(i6+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i6+x2.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x2, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x5, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x5, x2, x6, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x5.size(), x2.size(), x6.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x5, x2, x6, x4, x3, x1, x0);
}

void Task25::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x5 = b(3);
  const Index x6 = b(4);
  const Index x2 = b(5);
  const Index x3 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma115
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x3, x2, x6, x5, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x3, x2, x6, x5, x4, x1, x0), 0.0);
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i3+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i4+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))))];
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
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x2 == x4 && x5 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i3+x3.size()*(i4+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                  += (-2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x6);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i4+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x2, x6, x5, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4)))))];
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
    if (x5 == x6 && x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i3+x3.size()*(i3+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i3+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x2, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i3+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i3+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i4+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x5, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i4+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i6+x6.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x3, x2, x6, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x3.size(), x2.size(), x6.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x3, x2, x6, x5, x4, x1, x0);
}

void Task26::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma119
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x5, x4, x2, x3, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x4, x2, x3, x1, x0), 0.0);
  {
    if (x2 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x6) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x5, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x0, x5, x4, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i6+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3)))))];
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
    if (x2 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x2 == x4 && x5 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))))]
                  += (2.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x5 == x6 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x5, x0, x2, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i4+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i3)))))];
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
    if (x5 == x6 && x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (2.0) * i0data[i7+x7.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    if (x5 == x6) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x4, x2, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(3)->get_block(x7, x6, x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x6, x5, x4, x2, x3, x1, x0);
}

void Task27::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x3 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma122
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x3, x6, x5, x2, x4, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x3, x6, x5, x2, x4, x1, x0), 0.0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i4)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x6, x5, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4)))))];
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
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x6, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))];
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
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x3, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                    += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x6, x0, x2, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i0+x0.size()*(i2+x2.size()*(i4)))))];
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
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x4, x6, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x6, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x6, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x6, x5, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x3.size(), x6.size(), x5.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x7, x3, x6, x5, x2, x4, x1, x0);
}

void Task28::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x2 = b(3);
  const Index x6 = b(4);
  const Index x7 = b(5);
  const Index x3 = b(6);
  const Index x4 = b(7);
  // tensor label: Gamma550
  std::unique_ptr<double[]> odata(new double[out()->get_size(x7, x6, x2, x5, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x2, x5, x1, x0), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(4)->get_block(x4, x3);
  if (x2 == x6 && x7 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i5+x1.size()*(i0)))))]
                += (-1.0) * i0data[i4+x4.size()*(i0)] * fdata[i4+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x6, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i5+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i6+x6.size()*(i2+x2.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x2 == x5 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i3+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i6+x1.size()*(i0)))))]
                += (2.0) * i0data[i4+x4.size()*(i0)] * fdata[i4+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x0, x2, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i3+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i6+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i5)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x6, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (2.0) * i0data[i4+x4.size()*(i6+x6.size()*(i1+x1.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x2 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i3+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x6, x2, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i3+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i4+x4.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i4+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x5 && x4 == x6 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i3+x1.size()*(i0)))))]
                += (-2.0) * i0data[i7+x7.size()*(i0)] * fdata[i6+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x6 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i3+x1.size()*(i0)))))]
                += (1.0) * i0data[i7+x7.size()*(i0)] * fdata[i5+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i0)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x3 && x4 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i5)))] * fdata[i6+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x6 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i5+x1.size()*(i0)))))]
                += (1.0) * i0data[i7+x7.size()*(i0)] * fdata[i6+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x4 == x6 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i5+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i0)))] * fdata[i6+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i6+x1.size()*(i0)))))]
                += (-2.0) * i0data[i7+x7.size()*(i0)] * fdata[i5+x4.size()*(i3)];
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x1 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x0, x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i6+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i3)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x6, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x6 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i1+x1.size()*(i0)))] * fdata[i6+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x6 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i6+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x7, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x6, x2, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i5+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x6) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x7, x3, x2, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x7, x6, x2, x5, x1, x0);
}

void Task29::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x7 = b(2);
  const Index x2 = b(3);
  const Index x8 = b(4);
  const Index x9 = b(5);
  const Index x3 = b(6);
  const Index x4 = b(7);
  const Index x5 = b(8);
  const Index x6 = b(9);
  // tensor label: Gamma551
  std::unique_ptr<double[]> odata(new double[out()->get_size(x9, x8, x2, x7, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x9, x8, x2, x7, x1, x0), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(4)->get_block(x6, x5, x4, x3);
  if (x9 == x5 && x2 == x8 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i5+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x8, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i5+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x5 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x8);
    for (int i8 = 0; i8 != x8.size(); ++i8) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i8)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x8, x4, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i5+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x7 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i5+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (2.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x3 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i5+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i7)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x5 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x7, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i3+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x0, x4, x3, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i5+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x8, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i5+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (2.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x7, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i5+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x8, x4, x0, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i5+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x7 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x8);
    for (int i8 = 0; i8 != x8.size(); ++i8) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i3+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (2.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i8)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x4, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i3+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i7)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x0, x4, x8, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i8+x8.size()*(i2+x2.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x8, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i5+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (2.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x7, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i5+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x8, x4, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i5+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x7, x4, x8, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i3+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i8+x8.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x6, x8, x4, x3, x2, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i5 = 0; i5 != x5.size(); ++i5) {
                      odata[i5+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x8 && x1 == x7 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i6+x6.size()*(i0)] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x7 && x2 == x8 && x9 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i3+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x8, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i2+x2.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x5, x4, x8, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i8+x8.size()*(i2+x2.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3 && x2 == x7 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                odata[i3+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                  += (2.0) * i0data[i6+x6.size()*(i0)] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x7 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i3+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (2.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x0, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i2+x2.size()*(i7)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x5, x4, x0, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x8, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (2.0) * i0data[i6+x6.size()*(i8+x8.size()*(i1+x1.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x5, x4, x8, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i3+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (2.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i8+x8.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x6, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i3+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (-1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i1+x1.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x5, x4, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i3+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x9 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x6, x8, x2, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-1.0) * i0data[i6+x6.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x6, x5, x4, x8, x2, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i3 = 0; i3 != x3.size(); ++i3) {
                      odata[i3+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (-1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x3 && x4 == x8 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i5)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x8 && x4 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i5)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x7 && x4 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x5 && x2 == x8 && x4 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x2 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i5)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x1 == x7 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x7 && x4 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x6, x5, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x3 && x1 == x8 && x4 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i5)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x1 == x8 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i3)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x0, x6, x5, x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i3)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x6, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-2.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x6, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x3 && x2 == x5 && x6 == x8 && x4 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i9+x9.size()*(i0)] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x2 == x5 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x5 && x6 == x7 && x1 == x3 && x4 == x8) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i9+x9.size()*(i0)] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x5 && x4 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x6, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i6+x6.size()*(i7)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x1 == x3 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i8+x8.size()*(i4+x4.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x5 && x1 == x3 && x6 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x7, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i7+x7.size()*(i4+x4.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x3 && x4 == x5 && x2 == x7 && x6 == x8) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i9+x9.size()*(i0)] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x2 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i9+x9.size()*(i0)] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x2, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i7)))] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i2+x2.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x6, x5, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x2, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i5)))] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x0, x2, x7, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i7+x7.size()*(i6+x6.size()*(i5)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x2, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x1 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x2, x7, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i3+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i2+x2.size()*(i7+x7.size()*(i4+x4.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i9+x9.size()*(i0)] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                  += (1.0) * i0data[i9+x9.size()*(i0)] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x7, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i7+x7.size()*(i6+x6.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i4+x4.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x3 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i7)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x7 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i3)))] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x6, x0, x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i0+x0.size()*(i2+x2.size()*(i3)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x2, x7, x6, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i7+x7.size()*(i6+x6.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x4, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x1 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x0, x2, x7, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i5+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i7+x7.size()*(i4+x4.size()*(i3)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5 && x2 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                  += (1.0) * i0data[i9+x9.size()*(i0)] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x3 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x5 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x1 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x4, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i7+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x2 == x3 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                  += (-2.0) * i0data[i9+x9.size()*(i0)] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x3 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i5)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x5 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i0+x0.size()*(i4+x4.size()*(i3)))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x0, x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i3)))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x1 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x0, x2, x5, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i8+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x8, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i8 = 0; i8 != x8.size(); ++i8) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i8+x8.size()*(i1+x1.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i7+x7.size()*(i1+x1.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i5+x5.size()*(i1+x1.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x6, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-2.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i1+x1.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x7, x6, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x4, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x4, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i8+x8.size()*(i3+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x6, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x6, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i6+x6.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-2.0) * i0data[i9+x9.size()*(i8+x8.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x7, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i5+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (-2.0) * i0data[i9+x9.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x2 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i8+x8.size()*(i7+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (-2.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x9, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                    += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i1+x1.size()*(i0)))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i8+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x8, x2, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i8 = 0; i8 != x8.size(); ++i8) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x2, x7, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8 && x4 == x7) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x5, x2, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))] * fdata[i8+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x7) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x9, x8, x6, x5, x2, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i7 = 0; i7 != x7.size(); ++i7) {
                      odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x4 == x8) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x9, x3, x2, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i8 = 0; i8 != x8.size(); ++i8) {
                    odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                      += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x8) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x9, x3, x2, x7, x6, x5, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i8 = 0; i8 != x8.size(); ++i8) {
                      odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (1.0) * i0data[i9+x9.size()*(i3+x3.size()*(i2+x2.size()*(i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i8+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x9, x8, x2, x5, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i8 = 0; i8 != x8.size(); ++i8) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i7 = 0; i7 != x7.size(); ++i7) {
                      odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (1.0) * i0data[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x8) {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x9, x5, x2, x7, x4, x3, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i8 = 0; i8 != x8.size(); ++i8) {
                      odata[i9+x9.size()*(i8+x8.size()*(i2+x2.size()*(i7+x7.size()*(i1+x1.size()*(i0)))))]
                        += (1.0) * i0data[i9+x9.size()*(i5+x5.size()*(i2+x2.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))] * fdata[i8+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x9, x8, x2, x7, x1, x0);
}

void Task30::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma10
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x0, x1), 0.0);
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i2+x2.size()*(i1+x0.size()*(i1)))]
              += (2.0) * i0data[i3+x3.size()*(i2)];
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i2+x2.size()*(i2+x0.size()*(i1)))]
              += (-1.0) * i0data[i3+x3.size()*(i1)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x3, x2, x0, x1);
}

void Task31::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: Gamma12
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1), 0.0);
  {
    // rdm0 non-merged case
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        odata[i1+x0.size()*(i1)]  += 2.0 * i0data[0];
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->add_block(odata, x0, x1);
}

void Task32::Task_local::compute() {
  const Index x2 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma18
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x1, x0, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x1, x0, x2), 0.0);
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i1+x1.size()*(i1+x0.size()*(i2)))]
              += (1.0) * i0data[i3+x3.size()*(i2)];
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i1+x1.size()*(i2+x0.size()*(i2)))]
              += (-2.0) * i0data[i3+x3.size()*(i1)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x1, x0, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->add_block(odata, x3, x1, x0, x2);
}

void Task33::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: Gamma197
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))]
              += (2.0) * i0data[i2+x2.size()*(i1)];
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i1+x0.size()*(i3+x3.size()*(i3+x2.size()*(i1)))]  += 2.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i1+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))]
              += (-1.0) * i0data[i2+x2.size()*(i3)];
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            odata[i0+x0.size()*(i3+x3.size()*(i3+x2.size()*(i1)))]
              += (-1.0) * i0data[i0+x0.size()*(i1)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x0, x3, x2, x1);
}

void Task34::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x5 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma126
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x5, x0, x4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x3, x5, x0, x4, x2, x1), 0.0);
  {
    if (x2 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i5+x2.size()*(i1)))))]
                += (2.0) * i0data[i3+x3.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                += (-1.0) * i0data[i3+x3.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x0 == x1 && x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i5+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]  += -4.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x5);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                += (2.0) * i0data[i3+x3.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x3 == x4 && x2 == x5 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i4+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i5+x2.size()*(i1)))))]  += 2.0 * i0data[0];
          }
        }
      }
    }
  }
  {
    if (x2 == x5 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i5+x2.size()*(i1)))))]
                += (-1.0) * i0data[i3+x3.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x5);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                += (-1.0) * i0data[i2+x2.size()*(i5)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                += (2.0) * i0data[i2+x2.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x2, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i1+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x3.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                += (-4.0) * i0data[i2+x2.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i4+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (2.0) * i0data[i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x3.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                += (2.0) * i0data[i2+x2.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i5+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                += (2.0) * i0data[i0+x0.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x3, x5, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i5+x2.size()*(i1)))))]
                += (-1.0) * i0data[i0+x0.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i3+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i5+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x5, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i0+x0.size()*(i5+x5.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x3.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (2.0) * i0data[i0+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(3)->get_block(x3, x5, x0, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x3.size(), x5.size(), x0.size(), x4.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x3, x5, x0, x4, x2, x1);
}

void Task35::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma132
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x0, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x0, x3, x2, x1), 0.0);
  {
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i1+x0.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                += (2.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i1+x0.size()*(i3+x3.size()*(i4+x2.size()*(i1)))))]
                += (-1.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x0.size()*(i3+x3.size()*(i4+x2.size()*(i1)))))]
                += (2.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x0.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                += (-1.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x0, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i4+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i0+x0.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x0, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x5, x4, x0, x3, x2, x1);
}

void Task36::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma137
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x2, x3, x0, x1), 0.0);
  {
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i1+x0.size()*(i1)))))]
                += (4.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i1+x0.size()*(i1)))))]
                += (-2.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i1+x0.size()*(i1)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i3+x0.size()*(i1)))))]
                += (1.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i3+x0.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i4+x0.size()*(i1)))))]
                += (-2.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x2, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i3+x3.size()*(i4+x0.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i2+x2.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i3+x3.size()*(i0+x0.size()*(i1)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i3+x3.size()*(i0+x0.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x2, x3, x0, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x5, x4, x2, x3, x0, x1);
}

void Task37::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma155
  std::unique_ptr<double[]> odata(new double[out()->get_size(x2, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x2, x3, x0, x1), 0.0);
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x2.size()*(i3+x3.size()*(i3+x0.size()*(i1)))]
              += (1.0) * i0data[i2+x2.size()*(i1)];
          }
        }
      }
    }
  }
  {
    // rdm0 non-merged case
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block();
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i3+x2.size()*(i3+x3.size()*(i1+x0.size()*(i1)))]  += 4.0 * i0data[0];
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x2.size()*(i3+x3.size()*(i1+x0.size()*(i1)))]
              += (-2.0) * i0data[i2+x2.size()*(i3)];
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x2.size()*(i3+x3.size()*(i0+x0.size()*(i1)))]
              += (-2.0) * i0data[i0+x0.size()*(i1)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x2, x3, x0, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), x3.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x2, x3, x0, x1);
}

void Task38::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x3 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma176
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x3, x0, x4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x3, x0, x4, x2, x1), 0.0);
  {
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i1+x0.size()*(i4+x4.size()*(i3+x2.size()*(i1)))))]
                += (1.0) * i0data[i5+x5.size()*(i4)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i1+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                += (-2.0) * i0data[i5+x5.size()*(i3)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i1+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i3+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                += (1.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i3+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i4+x0.size()*(i4+x4.size()*(i3+x2.size()*(i1)))))]
                += (-2.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x0.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x0, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i4+x4.size()*(i3+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i0+x0.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x3, x0, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x3.size(), x0.size(), x4.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x5, x3, x0, x4, x2, x1);
}

void Task39::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma178
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x1, x0, x4, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x1, x0, x4, x3, x2), 0.0);
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i1+x1.size()*(i4+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i1+x1.size()*(i1+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i1+x1.size()*(i1+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i1+x1.size()*(i2+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (-2.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i1+x1.size()*(i2+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x0, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i1+x1.size()*(i0+x0.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i0+x0.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x1, x0, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x1.size(), x0.size(), x4.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x5, x1, x0, x4, x3, x2);
}

void Task40::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma179
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x3, x2, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x2, x0, x1), 0.0);
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x1, x3, x2);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i4+x0.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i1+x0.size()*(i1)))))]
                += (2.0) * i0data[i5+x5.size()*(i2)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x2);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1+x0.size()*(i1)))))]
                  += (2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i2+x0.size()*(i1)))))]
                += (-1.0) * i0data[i5+x5.size()*(i1)];
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i2+x0.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x0, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i0+x0.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i0+x0.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x5, x4, x3, x2, x0, x1);
}

void Task41::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x2 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma196
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x2, x4, x3, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x2, x4, x3, x0, x1), 0.0);
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i1+x0.size()*(i1)))))]
                  += (2.0) * i0data[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x1, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i2+x0.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i3+x0.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x2, x4, x3, x0, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x0.size(), x1.size());
  }
  out()->add_block(odata, x5, x2, x4, x3, x0, x1);
}

void Task42::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma552
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata.get(), out()->get_size(x0, x3), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x2, x1);
  // rdm0 merged case
  if (x0 == x1 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        odata[i1+x0.size()*(i3)]  += -2.0 * i0data[0] * fdata[i3+x2.size()*(i1)];
      }
    }
  }
  if (x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i0+x0.size()*(i3)]
            += (1.0) * i0data[i0+x0.size()*(i1)] * fdata[i3+x2.size()*(i1)];
        }
      }
    }
  }
  out()->add_block(odata, x0, x3);
}

void Task43::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x4 = b(5);
  // tensor label: Gamma554
  std::unique_ptr<double[]> odata(new double[out()->get_size(x0, x5)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5), 0.0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(4)->get_block(x4, x3, x2, x1);
  if (x0 == x1 && x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i1+x0.size()*(i5)]
              += (-2.0) * i0data[i4+x4.size()*(i3)] * fdata[i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x2 == x5 && x0 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x4, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x0.size()*(i5)]
              += (1.0) * i0data[i4+x4.size()*(i1)] * fdata[i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x2 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x4, x3, x0, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i0+x0.size()*(i5)]
                += (1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i1)))] * fdata[i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i1)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x0 == x1 && x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          odata[i1+x0.size()*(i5)]  += -2.0 * i0data[0] * fdata[i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))];
        }
      }
    }
  }
  if (x0 == x1 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i1+x0.size()*(i5)]
              += (1.0) * i0data[i2+x2.size()*(i3)] * fdata[i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x0 == x3 && x4 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i3+x0.size()*(i5)]
              += (-2.0) * i0data[i2+x2.size()*(i1)] * fdata[i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x0, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i0+x0.size()*(i5)]
              += (1.0) * i0data[i0+x0.size()*(i1)] * fdata[i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x0, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i0+x0.size()*(i5)]
                += (1.0) * i0data[i0+x0.size()*(i3+x3.size()*(i2+x2.size()*(i1)))] * fdata[i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
            }
          }
        }
      }
    }
  }
  out()->add_block(odata, x0, x5);
}

void Task44::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma24
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x4, x1, x3, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x1, x3, x2, x0), 0.0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i4+x2.size()*(i0)))))]
                += (-2.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
                  += (-2.0) * i0data[i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x1, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i4+x2.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->add_block(odata, x5, x4, x1, x3, x2, x0);
}

void Task45::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  // tensor label: Gamma25
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x3, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x3, x2, x0), 0.0);
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))]
              += (-2.0) * i0data[i2+x2.size()*(i0)];
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            odata[i1+x1.size()*(i3+x3.size()*(i3+x2.size()*(i0)))]
              += (1.0) * i0data[i1+x1.size()*(i0)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->add_block(odata, x1, x3, x2, x0);
}

void Task46::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma27
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0, x1, x2), 0.0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i0+x0.size()*(i2+x1.size()*(i2)))]
              += (2.0) * i0data[i3+x3.size()*(i0)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x3, x0, x1, x2);
}

void Task47::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma29
  std::unique_ptr<double[]> odata(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))]
              += (1.0) * i0data[i3+x3.size()*(i0)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, x3, x2, x1, x0);
}

void Task48::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma31
  std::unique_ptr<double[]> odata(new double[out()->get_size(x5, x0, x4, x3, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x3, x1, x2), 0.0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))]
                  += (2.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x5, x0, x4, x3, x1, x2);
}

void Task49::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma32
  std::unique_ptr<double[]> odata(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x0), 0.0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->add_block(odata, x1, x0);
}

#endif
