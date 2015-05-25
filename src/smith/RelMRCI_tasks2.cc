//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_tasks2.cc
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

#include <src/smith/RelMRCI_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task50::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma33
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x4, x3, x1, x2);
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x1, x2);
}

void Task51::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma34
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task52::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma219
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x5, x2, x4, x3, x0);
  {
    if (x2 == x5 && x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x1.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                += (-1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                += (1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i4+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x3.size()*(i0)))))]
                += (1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x2, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i4+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x3.size()*(i0)))))]
                += (-1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x2, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i5+x3.size()*(i0)))))]
                += (-1.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                  += (-1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i4+x3.size()*(i0)))))]
                += (1.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i5+x2.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                  += (1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i4+x3.size()*(i0)))))]
                  += (1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x2, x4, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i5+x3.size()*(i0)))))]
                  += (1.0) * i0data[i2+x2.size()*(i4+x4.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x5, x2, x4, x3, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x1.size(), x5.size(), x2.size(), x4.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, x1, x5, x2, x4, x3, x0);
}

void Task53::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma220
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x1, x5, x4, x0, x3, x2);
  {
    if (x3 == x6 && x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i6+x3.size()*(i2)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i0)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x4 == x6 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))))]
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
    if (x4 == x5 && x3 == x6 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i6+x3.size()*(i2)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x5, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i6+x3.size()*(i2)))))))]
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
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i5)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x3, x5, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x5, x4, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))];
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
    if (x4 == x6 && x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i1+x1.size()*(i2)))))];
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
    if (x4 == x5 && x3 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x1, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i6+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x2, x1, x5, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i6+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x3, x2, x1, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x1, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x7, x6, x1, x5, x4, x0, x3, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x1.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x7, x6, x1, x5, x4, x0, x3, x2);
}

void Task54::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma221
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x1, x5, x4, x3, x2, x0);
  {
    if (x2 == x6 && x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))))]
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
    if (x2 == x5 && x1 == x3 && x4 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x5, x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x1, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i5)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x1, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x1, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i6+x2.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i3)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x1, x5, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i6+x2.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x1, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x3, x1, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i6+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i1+x1.size()*(i5+x5.size()*(i2+x2.size()*(i0)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x7, x6, x1, x5, x4, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x1.size(), x5.size(), x4.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x7, x6, x1, x5, x4, x3, x2, x0);
}

void Task55::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma224
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x5, x4, x0, x3, x2);
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))]
                += (-1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x3.size()*(i2)))))]
                  += (1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i2+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))]
                += (1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x5, x4, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i2+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2, x1, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x5, x4, x0, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x1.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x1, x5, x4, x0, x3, x2);
}

void Task56::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma226
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x1, x5, x4, x3, x2, x0);
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))]
                += (1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                += (-1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i3+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
                += (-1.0) * i0data[i2+x2.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x5, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                += (1.0) * i0data[i1+x1.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                  += (1.0) * i0data[i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[i1+x1.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
                  += (1.0) * i0data[i1+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x5, x4, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x1.size(), x5.size(), x4.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x1, x5, x4, x3, x2, x0);
}

void Task57::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: Gamma225
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x4, x5, x3, x2, x1, x0);
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                += (-1.0) * i0data[i4+x4.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                += (1.0) * i0data[i3+x3.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x4, x5, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x4, x5, x3, x2, x1, x0);
}

void Task58::Task_local::compute() {
  const Index x6 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x7 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  const Index x4 = b(6);
  const Index x5 = b(7);
  // tensor label: Gamma234
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x1, x6);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(4)->get_block(x5, x4, x3, x2);
  if (x3 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i0+x0.size()*(i2+x1.size()*(i6)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i0+x0.size()*(i4+x1.size()*(i6)))]
                  += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x7.size()*(i0+x0.size()*(i6+x1.size()*(i6)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x5, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x1.size()*(i6)))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x5 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i7+x7.size()*(i0+x0.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i7+x7.size()*(i0)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i0+x0.size()*(i2+x1.size()*(i6)))]
                  += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i0+x0.size()*(i2+x1.size()*(i6)))]
                  += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i4)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x3, x6, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i7+x7.size()*(i0+x0.size()*(i2+x1.size()*(i6)))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i0+x0.size()*(i4+x1.size()*(i6)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x5, x6, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i0+x0.size()*(i4+x1.size()*(i6)))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i3+x3.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6)))]
                  += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x1, x6, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x1, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6)))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x7, x0, x1, x6, x5, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6)))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x7, x0, x1, x6);
}

void Task59::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma235
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x1, x4, x3, x2);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i0+x0.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x0, x1, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x1.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x1, x4, x3, x2);
}

void Task60::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma237
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x2, x4, x1, x3);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i0+x0.size()*(i4+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
                += (-1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x2, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i4+x4.size()*(i3+x1.size()*(i3)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i0+x0.size()*(i3+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x2, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i4+x4.size()*(i4+x1.size()*(i3)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x4);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i3+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x2.size()*(i4+x4.size()*(i1+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x0, x2, x4, x1, x3);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x2.size(), x4.size(), x1.size(), x3.size());
  }
  out()->put_block(odata, x5, x0, x2, x4, x1, x3);
}

void Task61::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma239
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x3, x4, x1, x2);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i0+x0.size()*(i4+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4+x4.size()*(i2+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x0, x3, x4, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x3.size(), x4.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x3, x4, x1, x2);
}

void Task62::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma238
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x3, x1, x4, x2, x0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i3+x1.size()*(i4+x4.size()*(i4+x2.size()*(i0)))))]
                += (-1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i3+x1.size()*(i4+x4.size()*(i2+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i3+x3.size()*(i4+x1.size()*(i4+x4.size()*(i3+x2.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x3, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x1.size()*(i4+x4.size()*(i2+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i4+x4.size()*(i3+x2.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i4+x4.size()*(i4+x2.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x3, x1, x4, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x3.size(), x1.size(), x4.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x5, x3, x1, x4, x2, x0);
}

void Task63::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma240
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x2, x1, x4, x3, x0);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x2, x1, x4, x3, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x2.size(), x1.size(), x4.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x1, x4, x3, x0);
}

void Task64::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x6 = b(2);
  const Index x7 = b(3);
  const Index x2 = b(4);
  const Index x3 = b(5);
  const Index x4 = b(6);
  const Index x5 = b(7);
  // tensor label: Gamma245
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x1, x0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(4)->get_block(x5, x4, x3, x2);
  if (x3 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x5, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i5+x5.size()*(i0)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x5, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x2, x5, x4, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x5 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i0)))]
                += (1.0) * i0data[i7+x7.size()*(i0)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x5, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0)))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x3, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i0)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x5, x4, x3, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  odata[i7+x7.size()*(i6+x6.size()*(i2+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i6+x6.size()*(i4+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x5, x0, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))]
                  += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i0)))] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x5, x2, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))] * fdata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x4, x3, x2, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))]
                    += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x7, x6, x5, x4, x3, x2, x1, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i1+x1.size()*(i0)))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))] * fdata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x7, x6, x1, x0);
}

void Task65::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma246
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x4, x3, x2, x1, x0);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i4+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x3, x2, x1, x0);
}

void Task66::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma250
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x4, x3, x0, x1, x2);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x1.size()*(i2)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2+x1.size()*(i2)))))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i1+x1.size()*(i2)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i1+x1.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x1, x2);
}

void Task67::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma256
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x6, x5, x1, x4, x3, x2);
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
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
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i4)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x3, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i4)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x4, x3, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i4+x3.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x2, x1, x4);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i5+x3.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2+x2.size()*(i1+x1.size()*(i4)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x6, x5, x1, x4, x3, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x1, x4, x3, x2);
}

void Task68::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma257
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x6, x5, x4, x3, x1, x2);
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x2, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2+x2.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x3);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i2)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x3, x1, x2);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i1+x1.size()*(i2)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x6, x5, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x4, x3, x1, x2);
}

void Task69::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x2 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma258
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x2, x4, x3, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x2, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x4, x3, x1, x0);
}

void Task70::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma282
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x4, x2, x1, x3);
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x3);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2+x2.size()*(i2+x1.size()*(i3)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x2);
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2+x2.size()*(i3+x1.size()*(i3)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x4, x2, x1, x3);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x2.size(), x1.size(), x3.size());
  }
  out()->put_block(odata, x5, x0, x4, x2, x1, x3);
}

void Task71::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma284
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x4, x5, x3, x2, x1, x0);
  {
    if (x3 == x6 && x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i6+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i2)))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x6 && x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
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
    if (x4 == x5 && x3 == x6 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i6+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                  += (1.0) * i0data[i7+x7.size()*(i0)];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x6 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
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
    if (x4 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i6+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x4, x5, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i6+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
    if (x4 == x6 && x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x4, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))];
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
    if (x4 == x5 && x3 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i6+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x2, x4, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i5+x5.size()*(i6+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i4+x4.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x3, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x5, x3, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x4.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x7, x6, x4, x5, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), x4.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x7, x6, x4, x5, x3, x2, x1, x0);
}

void Task72::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma303
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x2, x3, x4, x1, x0);
  {
    if (x1 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i4+x4.size()*(i4+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i2+x2.size()*(i4+x3.size()*(i4+x4.size()*(i2+x1.size()*(i0)))))]
                += (1.0) * i0data[i5+x5.size()*(i0)];
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x4);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i3+x3.size()*(i4+x4.size()*(i2+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i2+x2.size()*(i4+x3.size()*(i4+x4.size()*(i1+x1.size()*(i0)))))]
                  += (1.0) * i0data[i5+x5.size()*(i2+x2.size()*(i1+x1.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x2, x3, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x2.size(), x3.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x2, x3, x4, x1, x0);
}

void Task73::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x4 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma320
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x4, x6, x5, x3, x2, x1, x0);
  {
    if (x3 == x4 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x5);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i4+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
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
    if (x3 == x5 && x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x4, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
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
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x6, x5, x3, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i4+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i2)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i4+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x2, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i4+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x6, x0, x3, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i5+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i0+x0.size()*(i3+x3.size()*(i2)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x6, x5, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i4+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x6, x2, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i5+x3.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i2+x2.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x4, x6, x5, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x7, x4, x6, x5, x3, x2, x1, x0);
}

void Task74::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x2 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma321
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x2, x6, x5, x4, x3, x1, x0);
  {
    if (x1 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x6, x0, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x3);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x2, x6, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x6, x5, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x2, x6, x3, x1, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i2 = 0; i2 != x2.size(); ++i2) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i2+x2.size()*(i6+x6.size()*(i3+x3.size()*(i1+x1.size()*(i0)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x2, x6, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x2.size(), x6.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x7, x2, x6, x5, x4, x3, x1, x0);
}

void Task75::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma346
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x3, x4, x2, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x3, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i2+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x2);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i3+x1.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x3, x4, x2, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x3.size(), x4.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x3, x4, x2, x1, x0);
}

void Task76::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma52
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x3, x4, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x3.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x3, x4, x2, x1);
}

void Task77::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma53
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}

void Task78::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma54
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x2, x1);
}

void Task79::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma55
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task80::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x3 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma347
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x3, x5, x4, x0, x2, x1);
  {
    if (x3 == x6 && x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i4+x4.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i6+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i4+x4.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i3+x3.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i6+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x1, x3, x5, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i6+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
    if (x4 == x6 && x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x4, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
    if (x4 == x5 && x3 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x5, x4, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x3.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x3, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x5, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i6+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x3, x5, x4, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x3.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x6, x3, x5, x4, x0, x2, x1);
}

void Task81::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: Gamma348
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x4, x5, x3, x0, x2, x1);
  {
    if (x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x4.size()*(i5+x5.size()*(i5+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x4 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x4.size()*(i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))]
                  += (1.0) * i0data[i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x5, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x4.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x4, x5, x3, x0, x2, x1);
}

void Task82::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x6 = b(4);
  const Index x5 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma349
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x5, x6, x4, x3, x2, x1);
  {
    if (x2 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x1, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i6+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i6+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x6, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i4+x4.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x5.size()*(i6+x6.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x5, x6, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x0.size(), x5.size(), x6.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x5, x6, x4, x3, x2, x1);
}

void Task83::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x5 = b(3);
  const Index x6 = b(4);
  const Index x3 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma350
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x3, x6, x5, x4, x2, x1);
  {
    if (x3 == x6 && x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x4, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x6, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i4+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x6, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i4+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i6+x6.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x3, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i3+x3.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x3, x6, x5, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x0.size(), x3.size(), x6.size(), x5.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x3, x6, x5, x4, x2, x1);
}

void Task84::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma353
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x6) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x1, x5, x0, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i6+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3)))))];
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
    if (x4 == x6 && x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x5, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i6+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i5+x5.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i6+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x5, x0, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x4, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x6, x5, x0, x4, x3, x2, x1);
}

void Task85::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  const Index x6 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma354
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x6, x5, x4, x3, x0, x2, x1);
  {
    if (x3 == x6 && x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x5, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i6+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i5+x5.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x1, x5, x4, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i6+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x5, x4, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i5+x5.size()*(i4+x4.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i3+x3.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x5, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i1+x1.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x6, x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i6+x6.size()*(i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x3, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i6 = 0; i6 != x6.size(); ++i6) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i6+x6.size()*(i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x6, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x6.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x6, x5, x4, x3, x0, x2, x1);
}

void Task86::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma357
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x6, x5, x4, x3, x2, x1);
  {
    if (x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x1, x4, x3);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i5+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i3)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i0 = 0; i0 != x0.size(); ++i0) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x3, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x0, x6, x5, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x5, x4, x3, x2, x1);
}

void Task87::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x4 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma358
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x4, x6, x5, x3, x0, x2, x1);
  {
    if (x3 == x5 && x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x6, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i5+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i6+x6.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x1, x6, x5, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i4+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i4+x3.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                    += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x6, x1, x3, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i3+x3.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i1+x1.size()*(i3+x3.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i4+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x4, x6, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i4 = 0; i4 != x4.size(); ++i4) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i5+x5.size()*(i5+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i4+x4.size()*(i6+x6.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x4, x6, x5, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x4, x6, x5, x3, x0, x2, x1);
}

void Task88::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x5 = b(4);
  const Index x6 = b(5);
  const Index x3 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma359
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x3, x6, x5, x4, x0, x2, x1);
  {
    if (x2 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x6, x1, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i5+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x1, x6, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i3+x2.size()*(i1)))))))]
                    += (-1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i6+x6.size()*(i0)))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x1, x6, x5, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i3+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x3, x6, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))]
                      += (-1.0) * i0data[i7+x7.size()*(i3+x3.size()*(i6+x6.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x7, x3, x6, x5, x4, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,-1,1>(i0data, odata, x7.size(), x3.size(), x6.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x3, x6, x5, x4, x0, x2, x1);
}

void Task89::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma367
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x3, x4, x0, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x4, x0);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i3+x3.size()*(i4+x4.size()*(i0+x0.size()*(i3+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x3, x4, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x3, x4, x0, x2, x1);
}

void Task90::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x4 = b(2);
  const Index x5 = b(3);
  const Index x3 = b(4);
  const Index x6 = b(5);
  const Index x0 = b(6);
  const Index x7 = b(7);
  // tensor label: Gamma374
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x6, x3, x5, x4, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1, x5, x4);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i3+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1+x1.size()*(i5+x5.size()*(i4)))))];
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x3, x5, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i4+x4.size()*(i4+x2.size()*(i1)))))))]
                      += (1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i5+x5.size()*(i1)))))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x3, x5, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x3.size(), x5.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x7, x0, x6, x3, x5, x4, x2, x1);
}

void Task91::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x5 = b(2);
  const Index x6 = b(3);
  const Index x0 = b(4);
  const Index x7 = b(5);
  const Index x3 = b(6);
  const Index x4 = b(7);
  // tensor label: Gamma564
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x7, x0, x6, x5, x2, x1);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(3)->get_block(x4, x3);
  if (x6 == x3 && x2 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x1, x7, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i7+x7.size()*(i0+x0.size()*(i3+x6.size()*(i5+x5.size()*(i5+x2.size()*(i1)))))]
                  += (1.0) * i0data[i4+x4.size()*(i1+x1.size()*(i7+x7.size()*(i0)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x5, x7, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i7+x7.size()*(i0+x0.size()*(i3+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))]
                    += (1.0) * i0data[i4+x4.size()*(i5+x5.size()*(i7+x7.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))] * fdata[i4+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3 && x2 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i3+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i5+x2.size()*(i1)))))]
                  += (1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i6+x6.size()*(i1)))] * fdata[i4+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x7 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x0, x6, x5, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  odata[i3+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))]
                    += (1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))] * fdata[i4+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i3+x2.size()*(i1)))))]
                  += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i1)))] * fdata[i5+x4.size()*(i3)];
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))]
                    += (-1.0) * i0data[i7+x7.size()*(i0+x0.size()*(i6+x6.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))] * fdata[i5+x4.size()*(i3)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x7, x0, x6, x5, x2, x1);
}

void Task92::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x7 = b(2);
  const Index x8 = b(3);
  const Index x0 = b(4);
  const Index x9 = b(5);
  const Index x3 = b(6);
  const Index x4 = b(7);
  const Index x5 = b(8);
  const Index x6 = b(9);
  // tensor label: Gamma565
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x9, x0, x8, x7, x2, x1);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(3)->get_block(x6, x5, x4, x3);
  if (x8 == x5 && x2 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x7, x4, x1, x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i0+x0.size()*(i5+x8.size()*(i7+x7.size()*(i3+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i1+x1.size()*(i9+x9.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x3 && x2 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x1, x4, x7, x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i5+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i7+x7.size()*(i9+x9.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x3 && x9 == x5 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x6, x0, x4, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i5+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                    += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i1)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x8 == x5 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x6, x1, x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x9.size()*(i0+x0.size()*(i5+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                    += (1.0) * i0data[i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x5 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x1, x4, x3, x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i0+x0.size()*(i5+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i9+x9.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x0, x4, x3, x8, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i5+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i8+x8.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x0, x4, x1, x8, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i5+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i3+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i8+x8.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x1, x4, x0, x8, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i3+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i5+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i1+x1.size()*(i4+x4.size()*(i0+x0.size()*(i8+x8.size()*(i7)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x5 && x8 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x0, x4, x7, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i5+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x8 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x7, x4, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i3+x9.size()*(i0+x0.size()*(i5+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x6, x7, x4, x3, x9, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i5 = 0; i5 != x5.size(); ++i5) {
                      odata[i9+x9.size()*(i0+x0.size()*(i5+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i4+x4.size()*(i3+x3.size()*(i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
  if (x9 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x6, x0, x4, x3, x8, x7, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i5 = 0; i5 != x5.size(); ++i5) {
                      odata[i5+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
  if (x8 == x3 && x4 == x5 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x6, x1, x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i9+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                    += (1.0) * i0data[i6+x6.size()*(i1+x1.size()*(i9+x9.size()*(i0)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x3 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x5, x4, x1, x9, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i9 = 0; i9 != x9.size(); ++i9) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i9+x9.size()*(i0)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x6, x0, x8, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  odata[i3+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                    += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i8+x8.size()*(i1)))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x9 == x3 && x2 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x5, x4, x0, x8, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i3+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i7+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i8+x8.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x8 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x7, x9, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i7+x7.size()*(i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x8 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x6, x5, x4, x7, x9, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i3 = 0; i3 != x3.size(); ++i3) {
                      odata[i9+x9.size()*(i0+x0.size()*(i3+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i7+x7.size()*(i9+x9.size()*(i0+x0.size()*(i2+x2.size()*(i1)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
  if (x4 == x5 && x9 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x6, x0, x8, x7, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i3+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                      += (1.0) * i0data[i6+x6.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x6, x5, x4, x0, x8, x7, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i6 = 0; i6 != x6.size(); ++i6) {
                    for (int i3 = 0; i3 != x3.size(); ++i3) {
                      odata[i3+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (1.0) * i0data[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
  if (x4 == x7 && x2 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x9, x0, x8, x1, x6, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i3+x2.size()*(i1)))))]
                      += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i1+x1.size()*(i6+x6.size()*(i5)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x5 && x4 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x9, x0, x8, x3, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  for (int i7 = 0; i7 != x7.size(); ++i7) {
                    odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i5+x2.size()*(i1)))))]
                      += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i3+x3.size()*(i6+x6.size()*(i1)))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x9, x0, x8, x3, x6, x5, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i7 = 0; i7 != x7.size(); ++i7) {
                      odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i3+x3.size()*(i6+x6.size()*(i5+x5.size()*(i2+x2.size()*(i1)))))))] * fdata[i6+x6.size()*(i5+x5.size()*(i7+x4.size()*(i3)))];
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
  if (x2 == x3 && x6 == x7 && x4 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x9, x0, x8, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                for (int i5 = 0; i5 != x5.size(); ++i5) {
                  odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i3+x2.size()*(i1)))))]
                    += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i1)))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x6 == x7 && x2 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x9, x0, x8, x5, x4, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i3+x2.size()*(i1)))))]
                      += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i5+x5.size()*(i4+x4.size()*(i1)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x9, x0, x8, x1, x4, x3);
    for (int i3 = 0; i3 != x3.size(); ++i3) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i5+x2.size()*(i1)))))]
                      += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i1+x1.size()*(i4+x4.size()*(i3)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x9, x0, x8, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i8 = 0; i8 != x8.size(); ++i8) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i9 = 0; i9 != x9.size(); ++i9) {
                for (int i7 = 0; i7 != x7.size(); ++i7) {
                  for (int i5 = 0; i5 != x5.size(); ++i5) {
                    odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                      += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))] * fdata[i7+x6.size()*(i5+x5.size()*(i5+x4.size()*(i3)))];
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x9, x0, x8, x5, x4, x3, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              for (int i8 = 0; i8 != x8.size(); ++i8) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i9 = 0; i9 != x9.size(); ++i9) {
                    for (int i7 = 0; i7 != x7.size(); ++i7) {
                      odata[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i7+x7.size()*(i2+x2.size()*(i1)))))]
                        += (-1.0) * i0data[i9+x9.size()*(i0+x0.size()*(i8+x8.size()*(i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))))))] * fdata[i7+x6.size()*(i5+x5.size()*(i4+x4.size()*(i3)))];
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
  out()->put_block(odata, x9, x0, x8, x7, x2, x1);
}

void Task93::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma479
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x1, x4, x3, x2, x0);
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x4, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i4+x4.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x1, x4, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x1.size(), x4.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x5, x1, x4, x3, x2, x0);
}

void Task94::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  // tensor label: Gamma511
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x3, x1, x2, x0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x1.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x3, x1, x2, x0);
}

void Task95::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma534
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x1, x3, x2);
}

void Task96::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma570
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x3, x0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(2)->get_block(x2, x1);
  if (x3 == x1) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          odata[i1+x3.size()*(i0)]
            += (1.0) * i0data[i2+x2.size()*(i0)] * fdata[i2+x2.size()*(i1)];
        }
      }
    }
  }
  out()->put_block(odata, x3, x0);
}

void Task97::Task_local::compute() {
  const Index x0 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x4 = b(5);
  // tensor label: Gamma572
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(3)->get_block(x4, x3, x2, x1);
  if (x5 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x0, x2, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              odata[i3+x5.size()*(i0)]
                += (1.0) * i0data[i4+x4.size()*(i0+x0.size()*(i2+x2.size()*(i1)))] * fdata[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
            }
          }
        }
      }
    }
  }
  if (x2 == x3 && x5 == x1) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x4, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            odata[i1+x5.size()*(i0)]
              += (1.0) * i0data[i4+x4.size()*(i0)] * fdata[i4+x4.size()*(i3+x3.size()*(i3+x2.size()*(i1)))];
          }
        }
      }
    }
  }
  if (x5 == x1) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x4, x3, x2, x0);
    for (int i0 = 0; i0 != x0.size(); ++i0) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              odata[i1+x5.size()*(i0)]
                += (1.0) * i0data[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i0)))] * fdata[i4+x4.size()*(i3+x3.size()*(i2+x2.size()*(i1)))];
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, x5, x0);
}

void Task98::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma539
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x0, x2, x4, x3, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i4+x2.size()*(i4+x4.size()*(i3+x3.size()*(i1)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x2, x1);
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i4+x4.size()*(i4+x3.size()*(i1)))))]
                  += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i2+x2.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x2, x4, x3, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x2.size(), x4.size(), x3.size(), x1.size());
  }
  out()->put_block(odata, x5, x0, x2, x4, x3, x1);
}

void Task99::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x5 = b(5);
  // tensor label: Gamma540
  std::unique_ptr<std::complex<double>[]> odata = out()->move_block(x5, x4, x3, x1, x2, x0);
  {
    if (x2 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x1);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i1 = 0; i1 != x1.size(); ++i1) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x3.size()*(i1+x1.size()*(i4+x2.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i3+x3.size()*(i1)))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x3 == x4) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i4+x3.size()*(i1+x1.size()*(i2+x2.size()*(i0)))))]
                  += (-1.0) * i0data[i5+x5.size()*(i1+x1.size()*(i2+x2.size()*(i0)))];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x1, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), x1.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x3, x1, x2, x0);
}

#endif
