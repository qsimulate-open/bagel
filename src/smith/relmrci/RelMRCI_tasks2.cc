//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks2.cc
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

#include <src/smith/relmrci/RelMRCI_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

void Task50::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma31
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x4, x3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x5, x4, x3, x0, x2, x1);
}

void Task51::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma32
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x4, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x3, x2, x1), 0.0);
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
  out()->add_block(odata, x5, x0, x4, x3, x2, x1);
}

void Task52::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: Gamma33
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0, x2, x1), 0.0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->add_block(odata, x3, x0, x2, x1);
}

void Task53::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x3 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x0 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma230
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x6, x3, x5, x4, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x3, x5, x4, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x6, x3, x5, x4, x0, x2, x1);
}

void Task54::Task_local::compute() {
  const Index x4 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma231
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x4, x5, x3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x4, x5, x3, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x4, x5, x3, x0, x2, x1);
}

void Task55::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x5 = b(2);
  const Index x0 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma232
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x6, x5, x0, x4, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x0, x4, x3, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x6, x5, x0, x4, x3, x2, x1);
}

void Task56::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x0 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma233
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x6, x5, x4, x3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x6, x5, x4, x3, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x6, x5, x4, x3, x0, x2, x1);
}

void Task57::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma236
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x0, x6, x5, x4, x3, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x5, x4, x3, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x0, x6, x5, x4, x3, x2, x1);
}

void Task58::Task_local::compute() {
  const Index x7 = b(0);
  const Index x4 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x3 = b(4);
  const Index x0 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma237
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x4, x6, x5, x3, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x4, x6, x5, x3, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x4, x6, x5, x3, x0, x2, x1);
}

void Task59::Task_local::compute() {
  const Index x7 = b(0);
  const Index x3 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x0 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma238
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x3, x6, x5, x4, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x3, x6, x5, x4, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x3, x6, x5, x4, x0, x2, x1);
}

void Task60::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma245
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x3, x4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x3, x4, x2, x1), 0.0);
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
  out()->add_block(odata, x5, x0, x3, x4, x2, x1);
}

void Task61::Task_local::compute() {
  const Index x5 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma246
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x3, x4, x0, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x3, x4, x0, x2, x1), 0.0);
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
  out()->add_block(odata, x5, x3, x4, x0, x2, x1);
}

void Task62::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x3 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  const Index x2 = b(6);
  const Index x1 = b(7);
  // tensor label: Gamma253
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x0, x6, x3, x5, x4, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x3, x5, x4, x2, x1), 0.0);
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
  out()->add_block(odata, x7, x0, x6, x3, x5, x4, x2, x1);
}

void Task63::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  const Index x4 = b(6);
  const Index x3 = b(7);
  // tensor label: Gamma422
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x0, x6, x5, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x5, x2, x1), 0.0);
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
  if (x2 == x3 && x4 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i6 = 0; i6 != x6.size(); ++i6) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i3 = 0; i3 != x3.size(); ++i3) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
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
  out()->add_block(odata, x7, x0, x6, x5, x2, x1);
}

void Task64::Task_local::compute() {
  const Index x9 = b(0);
  const Index x0 = b(1);
  const Index x8 = b(2);
  const Index x7 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  const Index x6 = b(6);
  const Index x5 = b(7);
  const Index x4 = b(8);
  const Index x3 = b(9);
  // tensor label: Gamma423
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x9, x0, x8, x7, x2, x1)]);
  std::fill_n(odata.get(), out()->get_size(x9, x0, x8, x7, x2, x1), 0.0);
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
  if (x6 == x7 && x4 == x5 && x2 == x3) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x9, x0, x8, x1);
    for (int i1 = 0; i1 != x1.size(); ++i1) {
      for (int i8 = 0; i8 != x8.size(); ++i8) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i9 = 0; i9 != x9.size(); ++i9) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                for (int i3 = 0; i3 != x3.size(); ++i3) {
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
  out()->add_block(odata, x9, x0, x8, x7, x2, x1);
}

void Task65::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma317
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x4, x1, x3, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x1, x3, x2, x0), 0.0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i3 = 0; i3 != x3.size(); ++i3) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i5 = 0; i5 != x5.size(); ++i5) {
              odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i4+x2.size()*(i0)))))]
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
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              for (int i5 = 0; i5 != x5.size(); ++i5) {
                odata[i5+x5.size()*(i4+x4.size()*(i3+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))))]
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x3, x2, x0);
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x1, x0);
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
      std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x1, x3);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->add_block(odata, x5, x4, x1, x3, x2, x0);
}

void Task66::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: Gamma318
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x1, x3, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x1, x3, x2, x0), 0.0);
  {
    if (x1 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x2, x0);
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x1.size()*(i3+x3.size()*(i2+x2.size()*(i0)))]
              += (-1.0) * i0data[i2+x2.size()*(i0)];
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->add_block(odata, x1, x3, x2, x0);
}

void Task67::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma335
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x4, x3, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x3, x1, x2), 0.0);
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
  out()->add_block(odata, x5, x0, x4, x3, x1, x2);
}

void Task68::Task_local::compute() {
  const Index x5 = b(0);
  const Index x1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma336
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x1, x4, x3, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x1, x4, x3, x2, x0), 0.0);
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
  out()->add_block(odata, x5, x1, x4, x3, x2, x0);
}

void Task69::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: Gamma368
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x1, x2, x0), 0.0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x2, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x1.size(), x2.size(), x0.size());
  }
  out()->add_block(odata, x3, x1, x2, x0);
}

void Task70::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: Gamma363
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x0, x1, x2)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0, x1, x2), 0.0);
  {
    if (x1 == x2) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0);
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i0 = 0; i0 != x0.size(); ++i0) {
          for (int i3 = 0; i3 != x3.size(); ++i3) {
            odata[i3+x3.size()*(i0+x0.size()*(i2+x1.size()*(i2)))]
              += (1.0) * i0data[i3+x3.size()*(i0)];
          }
        }
      }
    }
  }
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->add_block(odata, x3, x0, x1, x2);
}

void Task71::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma391
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x4, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x1, x3, x2), 0.0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x5, x0, x4, x1, x3, x2);
}

void Task72::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: Gamma428
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x3, x0)]);
  std::fill_n(odata.get(), out()->get_size(x3, x0), 0.0);
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
  out()->add_block(odata, x3, x0);
}

void Task73::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma430
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0), 0.0);
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
  out()->add_block(odata, x5, x0);
}

void Task74::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  const Index x0 = b(5);
  // tensor label: Gamma396
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x4, x3, x1, x2, x0)]);
  std::fill_n(odata.get(), out()->get_size(x5, x4, x3, x1, x2, x0), 0.0);
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
  out()->add_block(odata, x5, x4, x3, x1, x2, x0);
}

void Task75::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x4 = b(4);
  const Index x1 = b(5);
  const Index x3 = b(6);
  const Index x2 = b(7);
  // tensor label: Gamma397
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x0, x6, x5, x4, x1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x5, x4, x1, x3, x2), 0.0);
  {
    if (x3 == x5) {
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x2, x4, x1);
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
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x7, x0, x6, x1, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x7, x0, x6, x5, x4, x1, x3, x2);
    sort_indices<0,1,2,3,4,5,6,7,1,1,1,1>(i0data, odata, x7.size(), x0.size(), x6.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size());
  }
  out()->add_block(odata, x7, x0, x6, x5, x4, x1, x3, x2);
}

void Task76::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  // tensor label: Gamma409
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x4, x2, x3, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x2, x3, x1), 0.0);
  {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x4, x2, x3, x1);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x2.size(), x3.size(), x1.size());
  }
  out()->add_block(odata, x5, x0, x4, x2, x3, x1);
}

void Task77::Task_local::compute() {
  const Index x0 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma412
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x5, x1, x4)]);
  std::fill_n(odata.get(), out()->get_size(x0, x5, x1, x4), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(4)->get_block(x3, x2);
  // rdm0 merged case
  if (x3 == x5 && x0 == x2 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          odata[i2+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4)))]  += -1.0 * i0data[0] * fdata[i5+x3.size()*(i2)];
        }
      }
    }
  }
  // rdm0 merged case
  if (x0 == x2 && x1 == x5 && x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          odata[i2+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4)))]  += 1.0 * i0data[0] * fdata[i4+x3.size()*(i2)];
        }
      }
    }
  }
  if (x3 == x4 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (-1.0) * i0data[i1+x1.size()*(i5)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x0 == x2 && x3 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i2+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (1.0) * i0data[i1+x1.size()*(i4)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x1 == x2 && x0 == x4 && x3 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          odata[i4+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]  += 1.0 * i0data[0] * fdata[i5+x3.size()*(i2)];
        }
      }
    }
  }
  if (x3 == x5 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
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
  if (x3 == x4 && x1 == x2 && x0 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i2 = 0; i2 != x2.size(); ++i2) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          odata[i5+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]  += -1.0 * i0data[0] * fdata[i4+x3.size()*(i2)];
        }
      }
    }
  }
  if (x0 == x5 && x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i5+x0.size()*(i5+x5.size()*(i1+x1.size()*(i4)))]
              += (1.0) * i0data[i1+x1.size()*(i2)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x5);
    for (int i5 = 0; i5 != x5.size(); ++i5) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]
              += (1.0) * i0data[i0+x0.size()*(i5)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x1 == x2 && x3 == x5) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            odata[i0+x0.size()*(i5+x5.size()*(i2+x1.size()*(i4)))]
              += (-1.0) * i0data[i0+x0.size()*(i4)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x5 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i0+x0.size()*(i5+x5.size()*(i4+x1.size()*(i4)))]
              += (1.0) * i0data[i0+x0.size()*(i2)] * fdata[i5+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x1 == x5 && x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i5 = 0; i5 != x5.size(); ++i5) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i0+x0.size()*(i5+x5.size()*(i5+x1.size()*(i4)))]
              += (-1.0) * i0data[i0+x0.size()*(i2)] * fdata[i4+x3.size()*(i2)];
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x5, x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x4, x0, x2);
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

void Task78::Task_local::compute() {
  const Index x0 = b(0);
  const Index x7 = b(1);
  const Index x1 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  const Index x3 = b(6);
  const Index x2 = b(7);
  // tensor label: Gamma413
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x7, x1, x6)]);
  std::fill_n(odata.get(), out()->get_size(x0, x7, x1, x6), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(5)->get_block(x5, x4, x3, x2);
  if (x0 == x2 && x1 == x6 && x3 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x2 && x3 == x6 && x1 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x4 && x1 == x6 && x3 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i4+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x0 == x4 && x1 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i4+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x6 && x1 == x2 && x3 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i6+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x1 == x4 && x3 == x7 && x0 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x6 && x3 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i4)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x1 == x4 && x0 == x7 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i7+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i2)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x7 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i1+x1.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                  += (1.0) * i0data[i5+x5.size()*(i4+x4.size()*(i0+x0.size()*(i2)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x1 == x7 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x4, x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i5 = 0; i5 != x5.size(); ++i5) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
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
  if (x0 == x2 && x1 == x4 && x5 == x7 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]  += 1.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x4 && x3 == x6 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i5+x5.size()*(i7)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x1 == x4 && x3 == x7 && x5 == x6 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]  += -1.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x4 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i6)] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x2 && x5 == x6 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x0 == x2 && x1 == x4 && x5 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i2+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x0 == x2 && x3 == x4 && x5 == x7 && x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]  += -1.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x6 && x5 == x7 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x0 == x2 && x3 == x4 && x1 == x7 && x5 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]  += 1.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x1 == x7 && x0 == x2 && x5 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i2+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x7);
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
  if (x3 == x4 && x5 == x7 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (1.0) * i0data[i1+x1.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x4);
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
  if (x0 == x2 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x7, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (1.0) * i0data[i1+x1.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x6, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i2+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x0 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x3, x7, x1, x4);
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
  if (x0 == x2 && x5 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x6, x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
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
  if (x1 == x2 && x0 == x4 && x5 == x7 && x3 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += -1.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x2 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i5+x5.size()*(i7)] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x6 && x3 == x7 && x1 == x2 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i7 = 0; i7 != x7.size(); ++i7) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += 1.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x2 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x6);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x7);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x6);
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i6)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x6 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x7 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i7+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x6 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (1.0) * i0data[i1+x1.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x0 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x7, x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x6, x5, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x7, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x6, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i4+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (1.0) * i0data[i1+x1.size()*(i6+x6.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  // rdm0 merged case
  if (x5 == x7 && x3 == x4 && x1 == x2 && x0 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            odata[i6+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += 1.0 * i0data[0] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x2 && x0 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x4);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i6+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x0 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x4, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block();
    for (int i6 = 0; i6 != x6.size(); ++i6) {
      for (int i4 = 0; i4 != x4.size(); ++i4) {
        for (int i2 = 0; i2 != x2.size(); ++i2) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]  += -1.0 * i0data[0] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x2 && x0 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i2 = 0; i2 != x2.size(); ++i2) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i3+x3.size()*(i4)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x1 == x4 && x0 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (-1.0) * i0data[i3+x3.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x0 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x1, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i1 = 0; i1 != x1.size(); ++i1) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                += (1.0) * i0data[i1+x1.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x0 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i1 = 0; i1 != x1.size(); ++i1) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i7 = 0; i7 != x7.size(); ++i7) {
                odata[i7+x0.size()*(i7+x7.size()*(i1+x1.size()*(i6)))]
                  += (1.0) * i0data[i1+x1.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x7);
    for (int i7 = 0; i7 != x7.size(); ++i7) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i0+x0.size()*(i7)] * fdata[i6+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x3 == x4 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x6);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i6 = 0; i6 != x6.size(); ++i6) {
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                += (1.0) * i0data[i0+x0.size()*(i4)] * fdata[i7+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x6 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x7, x5, x4);
    for (int i4 = 0; i4 != x4.size(); ++i4) {
      for (int i5 = 0; i5 != x5.size(); ++i5) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                odata[i0+x0.size()*(i7+x7.size()*(i2+x1.size()*(i6)))]
                  += (1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i5+x5.size()*(i4)))] * fdata[i5+x5.size()*(i4+x4.size()*(i6+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x7 && x1 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x4);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x6, x5, x4);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x7, x3, x4);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x3, x6, x0, x4);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x7, x5, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i6 = 0; i6 != x6.size(); ++i6) {
          for (int i7 = 0; i7 != x7.size(); ++i7) {
            for (int i4 = 0; i4 != x4.size(); ++i4) {
              odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                += (1.0) * i0data[i0+x0.size()*(i2)] * fdata[i6+x5.size()*(i4+x4.size()*(i7+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x3 == x7 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x5, x6, x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x7, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              for (int i4 = 0; i4 != x4.size(); ++i4) {
                odata[i0+x0.size()*(i7+x7.size()*(i4+x1.size()*(i6)))]
                  += (1.0) * i0data[i0+x0.size()*(i7+x7.size()*(i3+x3.size()*(i2)))] * fdata[i6+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x4) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x6, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i0 = 0; i0 != x0.size(); ++i0) {
        for (int i7 = 0; i7 != x7.size(); ++i7) {
          for (int i4 = 0; i4 != x4.size(); ++i4) {
            for (int i6 = 0; i6 != x6.size(); ++i6) {
              odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                += (1.0) * i0data[i0+x0.size()*(i2)] * fdata[i7+x5.size()*(i4+x4.size()*(i4+x3.size()*(i2)))];
            }
          }
        }
      }
    }
  }
  if (x5 == x7 && x1 == x6) {
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x4, x3, x2);
    for (int i2 = 0; i2 != x2.size(); ++i2) {
      for (int i3 = 0; i3 != x3.size(); ++i3) {
        for (int i4 = 0; i4 != x4.size(); ++i4) {
          for (int i0 = 0; i0 != x0.size(); ++i0) {
            for (int i7 = 0; i7 != x7.size(); ++i7) {
              for (int i6 = 0; i6 != x6.size(); ++i6) {
                odata[i0+x0.size()*(i7+x7.size()*(i6+x1.size()*(i6)))]
                  += (1.0) * i0data[i0+x0.size()*(i4+x4.size()*(i3+x3.size()*(i2)))] * fdata[i7+x5.size()*(i4+x4.size()*(i3+x3.size()*(i2)))];
              }
            }
          }
        }
      }
    }
  }
  if (x5 == x6 && x3 == x4 && x1 == x7) {
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x4, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x7, x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x6, x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x0, x4, x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x0, x7, x5, x4, x1, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(2)->get_block(x1, x4, x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x1, x6, x5, x4, x0, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x0, x7, x1, x4, x3, x2);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(3)->get_block(x1, x6, x0, x4, x3, x2);
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

void Task79::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // scalar
  // tensor label: Gamma424
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(1)->get_block(x1, x0);
  out()->add_block(odata);
}

void Task80::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // scalar
  // tensor label: Gamma426
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size()]);
  std::fill_n(odata.get(), out()->get_size(), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(2)->get_block(x3, x2, x1, x0);
  out()->add_block(odata);
}

void Task81::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  // tensor label: Gamma432
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x5, x0, x4, x1)]);
  std::fill_n(odata.get(), out()->get_size(x5, x0, x4, x1), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(2)->get_block(x3, x2);
  if (x4 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x1, x5, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x3, x0, x4, x1);
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

void Task82::Task_local::compute() {
  const Index x7 = b(0);
  const Index x0 = b(1);
  const Index x6 = b(2);
  const Index x1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  const Index x3 = b(6);
  const Index x2 = b(7);
  // tensor label: Gamma433
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x7, x0, x6, x1)]);
  std::fill_n(odata.get(), out()->get_size(x7, x0, x6, x1), 0.0);
  // associated with merged
  std::unique_ptr<std::complex<double>[]> fdata = in(3)->get_block(x5, x4, x3, x2);
  if (x7 == x4 && x6 == x2) {
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x3, x1);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x3, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x1, x3, x2, x7, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x0, x3, x2, x6, x1);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x1, x7, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x1, x7, x0);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x0, x6, x1);
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
    std::unique_ptr<std::complex<double>[]> i0data = in(1)->get_block(x5, x4, x3, x0, x6, x1);
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

void Task84::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x1, c1, x0)]);
  std::fill_n(odata.get(), out()->get_size(c2, x1, c1, x0), 0.0);
  {
    // tensor label: I0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, c1, x0, x1);
    sort_indices<0,3,1,2,1,1,1,1>(i0data, odata, c2.size(), c1.size(), x0.size(), x1.size());
  }
  {
    // tensor label: I0
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, c2, x1, x0);
    sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, c1.size(), c2.size(), x1.size(), x0.size());
  }
  out()->add_block(odata, c2, x1, c1, x0);
}

void Task85::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma0
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x3, x1, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      // tensor label: I1
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, c1, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, c1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c1.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task86::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x3, c3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x3, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,-2,1>(i1data, i1data_sorted, c2.size(), c3.size());
    zgemm3m_("T", "N", c1.size()*x3.size()*x2.size(), c2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task87::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c3, a4, c1, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c3, a4, c1, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x3.size());
      // tensor label: I214
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task88::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I214
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(a4, x2, c2, c3)]);
  std::fill_n(odata.get(), out()->get_size(a4, x2, c2, c3), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a4, x2, c2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), c2.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(c2, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c2.size(), x2.size(), a4.size(), c3.size());
  }
  out()->add_block(odata, a4, x2, c2, c3);
}

void Task89::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x2, a4, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x2.size(), a4.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), c2.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), c2.size(), x2.size());
  out()->add_block(odata, c2, c1, x3, x2);
}

void Task90::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I4
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x0, x1, x2, c1);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x0, x1, x2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), c1.size());
    zgemm3m_("T", "N", c2.size(), x0.size()*x1.size()*c1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), x1.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task91::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  // tensor label: I4
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: Gamma1
        std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(x5, x4, x0, x3, x1, x2);
        std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(x5, x4, x0, x3, x1, x2)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
        zgemm3m_("T", "N", c1.size(), x0.size()*x1.size()*x2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), x1.size(), x2.size());
  out()->add_block(odata, x0, x1, x2, c1);
}

void Task92::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x1, x3, x0, x2);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      // tensor label: I7
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, c2, x3, x2);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, c2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c2.size(), x3.size(), x2.size());
      zgemm3m_("T", "N", x1.size()*x0.size(), c1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task93::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I7
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(a3, x2);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c1, a3, c2, x3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c1, a3, c2, x3)]);
    sort_indices<1,0,2,3,0,1,-1,1>(i1data, i1data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    zgemm3m_("T", "N", x2.size(), c1.size()*c2.size()*x3.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), c1.size(), c2.size(), x3.size());
  out()->add_block(odata, c1, c2, x3, x2);
}

void Task94::Task_local::compute() {
  const Index c1 = b(0);
  const Index c2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I7
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c1, c2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c2, x3, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: v2
      std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(a4, x2, c2, c3);
      std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(a4, x2, c2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c2.size(), c3.size());
      zgemm3m_("T", "N", c1.size()*x3.size(), x2.size()*c2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<0,3,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), c2.size());
  out()->add_block(odata, c1, c2, x3, x2);
}

void Task95::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma58
          std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x5, x1, x4, x3, x2);
          std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x5, x1, x4, x3, x2)]);
          sort_indices<1,3,4,5,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x1.size(), x4.size(), x3.size(), x2.size());
          // tensor label: I183
          std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x3, x2, c1, x5, x4);
          std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x3, x2, c1, x5, x4)]);
          sort_indices<4,5,1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c1.size(), x5.size(), x4.size());
          zgemm3m_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task96::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I183
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x5, c3, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
    // tensor label: I184
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, c3, x3, x2);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, c3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size(), x3.size(), x2.size());
    zgemm3m_("T", "N", c1.size()*x5.size()*x4.size(), c2.size()*x3.size()*x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, c1, x5, x4);
}

void Task97::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I184
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c3, x3, x2)]);
  std::fill_n(odata.get(), out()->get_size(c2, c3, x3, x2), 0.0);
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c2, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(0)->get_block(x3, x2, c2, c3);
    sort_indices<2,3,0,1,1,1,-1,1>(i1data, odata, x3.size(), x2.size(), c2.size(), c3.size());
  }
  out()->add_block(odata, c2, c3, x3, x2);
}

void Task98::Task_local::compute() {
  const Index c2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I0
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c1, x0, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma59
          std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(x0, x5, x2, x4, x1, x3);
          std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(x0, x5, x2, x4, x1, x3)]);
          sort_indices<1,2,3,5,0,4,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size(), x2.size(), x4.size(), x1.size(), x3.size());
          // tensor label: I186
          std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x3, x2, c1, x5, x4);
          std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x3, x2, c1, x5, x4)]);
          sort_indices<4,2,5,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c1.size(), x5.size(), x4.size());
          zgemm3m_("T", "N", x0.size()*x1.size(), c2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x0.size()*x1.size());
        }
      }
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), c1.size());
  out()->add_block(odata, c2, c1, x0, x1);
}

void Task99::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I186
  std::unique_ptr<std::complex<double>[]> odata(new std::complex<double>[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  std::unique_ptr<std::complex<double>[]> odata_sorted(new std::complex<double>[out()->get_size(c2, x3, x2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x2, c1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<std::complex<double>[]> i0data = in(0)->get_block(c1, x5, c3, x4);
    std::unique_ptr<std::complex<double>[]> i0data_sorted(new std::complex<double>[in(0)->get_size(c1, x5, c3, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c3.size(), x4.size());
    // tensor label: v2
    std::unique_ptr<std::complex<double>[]> i1data = in(1)->get_block(c2, x3, x2, c3);
    std::unique_ptr<std::complex<double>[]> i1data_sorted(new std::complex<double>[in(1)->get_size(c2, x3, x2, c3)]);
    sort_indices<3,0,1,2,0,1,-1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x2.size(), c3.size());
    zgemm3m_("T", "N", c1.size()*x5.size()*x4.size(), c2.size()*x3.size()*x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), c2.size(), x3.size(), x2.size());
  out()->add_block(odata, c2, x3, x2, c1, x5, x4);
}

#endif
