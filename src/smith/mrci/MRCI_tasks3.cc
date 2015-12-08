//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks3.cc
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

#include <src/smith/MRCI_tasks3.h>

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
  // tensor label: Gamma531
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x5, x4, x1, x3, x2);
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
  out()->put_block(odata, x7, x0, x6, x5, x4, x1, x3, x2);
}

void Task101::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma543
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x2, x3, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x2, x3, x1);
}

void Task102::Task_local::compute() {
  const Index x4 = b(0);
  const Index x1 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma545
  std::unique_ptr<double[]> odata = out()->move_block(x0, x5, x1, x4);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(4)->get_block(x3, x2);
  // rdm0 merged case
  if (x3 == x5 && x0 == x2 && x1 == x4) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
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
  if (x3 == x5 && x0 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
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
  if (x3 == x4 && x0 == x5 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
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
  out()->put_block(odata, x0, x5, x1, x4);
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
  // tensor label: Gamma546
  std::unique_ptr<double[]> odata = out()->move_block(x0, x7, x1, x6);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(5)->get_block(x5, x4, x3, x2);
  if (x1 == x6 && x0 == x2 && x3 == x7) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (-4.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x1 == x7 && x0 == x2 && x3 == x6) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
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
  if (x0 == x2 && x5 == x6 && x1 == x4 && x3 == x7) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
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
  if (x5 == x6 && x1 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
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
  if (x5 == x6 && x3 == x7 && x0 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
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
  if (x5 == x6 && x3 == x7 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += 2.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x2 && x0 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(x5, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
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
  out()->put_block(odata, x0, x7, x1, x6);
}

void Task104::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // scalar
  // tensor label: Gamma555
  std::unique_ptr<double[]> odata = out()->move_block();
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  out()->put_block(odata);
}

void Task105::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // scalar
  // tensor label: Gamma557
  std::unique_ptr<double[]> odata = out()->move_block();
  // associated with merged
  std::unique_ptr<double[]> fdata = in(2)->get_block(x3, x2, x1, x0);
  out()->put_block(odata);
}

void Task106::Task_local::compute() {
  const Index x1 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  // tensor label: Gamma563
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
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
  out()->put_block(odata, x5, x0, x4, x1);
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
  // tensor label: Gamma564
  std::unique_ptr<double[]> odata = out()->move_block(x7, x0, x6, x1);
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
  out()->put_block(odata, x7, x0, x6, x1);
}

void Task109::compute_() {
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c2, c1, x0, x1") + (*ta1_)("c1, c2, x1, x0");
}

void Task110::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c2, c1, x3, x2");
}

void Task111::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("c1, x3, c3, x2") * (*ta2_)("c2, c3") * (-2);
}

void Task112::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("c3, a4, c1, x3") * (*ta2_)("a4, x2, c2, c3");
}

void Task113::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a4, x2, c2, c3") += (*ta1_)("a4, x2, c2, c3") * (-1)
     + (*ta1_)("c2, x2, a4, c3") * 2;
}

void Task114::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("c1, a4, c3, x3") * (*ta2_)("c2, x2, a4, c3") * (-1);
}

void Task115::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c2, x2") * (*ta2_)("c1, x0, x1, x2");
}

void Task116::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x0, x1, x2") += (*ta1_)("x5, x4, x0, x3, x1, x2") * (*ta2_)("x5, x4, c1, x3");
}

void Task117::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("x2, c1, c2, x3");
}

void Task118::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c1, c2, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, x2") * (-1);
}

void Task119::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x2, c1, c2, x3") += (*ta1_)("c1, a4, c3, x3") * (*ta2_)("a4, x2, c2, c3");
}

void Task120::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x5, x1, x4, x3, x2") * (*ta2_)("c2, x3, x2, c1, x5, x4");
}

void Task121::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, c1, x5, x4") += (*ta1_)("c1, x5, c3, x4") * (*ta2_)("c2, c3, x3, x2");
}

void Task122::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c3, x3, x2") += (*ta1_)("c2, c3, x3, x2") * (-1)
     + (*ta1_)("x3, x2, c2, c3") * (-1);
}

void Task123::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x5, x2, x4, x1, x3") * (*ta2_)("c2, x3, x2, c1, x5, x4");
}

void Task124::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, c1, x5, x4") += (*ta1_)("c1, x5, c3, x4") * (*ta2_)("c2, x3, x2, c3") * (-1);
}

void Task125::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x5, x3, x4, x1, x2") * (*ta2_)("x3, c2, x2, c1, x5, x4");
}

void Task126::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c2, x2, c1, x5, x4") += (*ta1_)("c1, x5, c3, x4") * (*ta2_)("x3, c3, c2, x2");
}

void Task127::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x7, x6, c1, x5") * (*ta2_)("c2, x7, x6, x0, x5, x1");
}

void Task128::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x6, x0, x5, x1") += (*ta1_)("x7, x6, x0, x5, x1, x4, x3, x2") * (*ta2_)("c2, x4, x3, x2") * 0.5;
}

void Task129::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, x7, x6, x0, x5, x1") += (*ta1_)("x7, x6, x0, x5, x4, x3, x1, x2") * (*ta2_)("x4, x3, c2, x2") * 0.5;
}

void Task130::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c1, x2, c2, c3") * (*ta2_)("c3, x1, x0, x2");
}

void Task131::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c3, x1, x0, x2") += (*ta1_)("x5, x4, x1, x3, x0, x2") * (*ta2_)("x5, x4, c3, x3");
}

void Task132::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c1, a3, c2, x5") * (*ta2_)("a3, x1, x5, x0");
}

void Task133::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x5, x0") += (*ta1_)("x1, x5, x0, x4, x3, x2") * (*ta2_)("a3, x4, x3, x2") * (-0.5);
}

void Task134::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("a3, x1, x5, x0") += (*ta1_)("x1, x5, x4, x3, x0, x2") * (*ta2_)("x4, x3, a3, x2") * (-0.5);
}

void Task135::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x5, x2, x0, x4, x1, x3") * (*ta2_)("c2, x3, x2, x5, c1, x4");
}

void Task136::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, x5, c1, x4") += (*ta1_)("x5, a3, c1, x4") * (*ta2_)("c2, x3, a3, x2");
}

void Task137::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x5, x4, x1, x3, x0, x2") * (*ta2_)("c2, x3, x2, c1, x5, x4");
}

void Task138::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c2, x3, x2, c1, x5, x4") += (*ta1_)("c1, a3, x5, x4") * (*ta2_)("c2, x3, a3, x2") * (-1);
}

void Task139::compute_() {
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
}

void Task140::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x3, x4, x1, x0") * (*ta2_)("x3, c1, x5, x4");
}

void Task141::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c1, x5, c2, x4") * (*ta2_)("x3, c2") * 2;
}

void Task142::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c2, x5, c3, x4") * (*ta2_)("x3, c3, c1, c2") * (-2);
}

void Task143::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c2, a3, c1, x5") * (*ta2_)("a3, x4, x3, c2") * 0.5;
}

void Task144::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3");
}

void Task145::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("c1, c2") * (-1);
}

void Task146::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("a2, x3");
}

void Task147::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, a3, c2, x4") * (*ta2_)("c1, x3, a3, c2") * (-1);
}

void Task148::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c2, a3, x5, x4") * (*ta2_)("a3, x3, c1, c2");
}

void Task149::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  (*ta0_)("a3, x3, c1, c2") += (*ta1_)("a3, x3, c1, c2") * (-1)
     + (*ta1_)("c1, x3, a3, c2") * 2;
}

#endif
