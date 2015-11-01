//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks2.cc
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

#include <src/smith/CASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task50::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma298
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x3, x2, x0);
  {
    if (x2 == x4 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x1.size()*(ix3+x3.size()*(ix4+x2.size()*(ix0))))))]
                  += (-2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))))]
                    += (-2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x3, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x1, x3);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix4+x2.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x3, x2, x0);
}

void Task51::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x0 = b(4);
  // tensor label: Gamma299
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x2, x0);
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x1.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))]
                += (-2.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix1+x1.size()*(ix3+x3.size()*(ix3+x2.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x3, x2, x0);
}

void Task52::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma301
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x1, x4);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x3, x2);
  if (x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x1.size()*(ix4))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix2))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x1.size()*(ix4))))]
                += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))] * fdata[ix4+x3.size()*(ix2)];
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x4);
    for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x1.size()*(ix4))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x1, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix2))))] * fdata[ix4+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x0, x1, x4, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix1+x1.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x0, x1, x4);
}

void Task53::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: Gamma302
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x1.size()*(ix2))))]
                += (2.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,-1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task54::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma304
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x3, x2);
  if (x1 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix2))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4 && x1 == x2) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x1.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0))] * fdata[ix4+x3.size()*(ix2)];
            }
          }
        }
      }
    }
  }
  if (x1 == x2) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x3, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix0))))] * fdata[ix3+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  if (x3 == x4) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x2, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))]
                  += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1+x1.size()*(ix0))))] * fdata[ix4+x3.size()*(ix2)];
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x3, x2, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x1.size()*(ix0))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1+x1.size()*(ix0))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task55::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: Gamma305
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix2+x1.size()*(ix0))))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0))];
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task56::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma307
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x1, x2);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x1.size()*(ix2))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x2);
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix2))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix2))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x1, x2);
}

void Task57::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma308
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task58::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma321
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x2, x4, x3, x1, x0);
  {
    if (x1 == x2) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x1.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x2, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix2+x2.size()*(ix4+x4.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x2, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x2, x4, x3, x1, x0);
}

void Task59::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma326
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x3, x4, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4+x4.size()*(ix4+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix3+x3.size()*(ix1))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x3.size()*(ix4+x4.size()*(ix2+x2.size()*(ix1))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x3, x4, x2, x1);
}

void Task60::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma327
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    if (x2 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1, x3, x0);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix0+x0.size()*(ix4+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix3+x3.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}

void Task61::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x5 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  const Index x4 = b(7);
  const Index x3 = b(8);
  // tensor label: Gamma328
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x7, x0, x6, x5, x2, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(3)->get_block(x4, x3);
  if (x2 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x1, x4, x3);
    for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix5+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix1+x1.size()*(ix4+x4.size()*(ix3))))))] * fdata[ix4+x4.size()*(ix3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5 && x2 == x3) {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x7, x0, x6, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix3+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix1))))] * fdata[ix5+x4.size()*(ix3)];
                }
              }
            }
          }
        }
      }
    }
  }
  if (x2 == x3) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x5, x4, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix3+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1))))))] * fdata[ix4+x4.size()*(ix3)];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (x4 == x5) {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x7, x0, x6, x3, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                  for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                    odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))]
                      += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))] * fdata[ix5+x4.size()*(ix3)];
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
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x7, x0, x6, x5, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
          for (int ix6 = 0; ix6 != x6.size(); ++ix6) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix7 = 0; ix7 != x7.size(); ++ix7) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix7+x7.size()*(ix0+x0.size()*(ix6+x6.size()*(ix5+x5.size()*(ix2+x2.size()*(ix1))))))];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x7, x0, x6, x5, x2, x1);
}

void Task62::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma329
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,1>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task63::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma330
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task64::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: Gamma339
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
          odata[ici0]
            += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))] * fdata[ix1+x1.size()*(ix0)];
        }
      }
    }
  }
  out()->put_block(odata, ci0);
}

void Task65::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: Gamma351
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
      for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
        for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
              odata[ici0+ci0.size()*(ix3+x3.size()*(ix0))]
                += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))] * fdata[ix2+x2.size()*(ix1)];
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task66::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x3 = b(5);
  const Index x2 = b(6);
  // tensor label: Gamma362
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  // associated with merged
  std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
  {
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1, x3, x2);
    for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
      for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))]
                    += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                }
              }
            }
          }
        }
      }
    }
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}

void Task67::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: Gamma377
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x5, x4, x3, x1, x0);
  {
    if (x2 == x5 && x1 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x1 == x5) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix5+x1.size()*(ix0))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x5) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x3, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix5+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix2+x2.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x3, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x1 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x2, x5, x4, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix3+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x4, x5, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix3+x2.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix4+x4.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x2, x3, x1, x0);
      for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
        for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix2+x2.size()*(ix5+x5.size()*(ix5+x4.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix2+x2.size()*(ix3+x3.size()*(ix1+x1.size()*(ix0))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x5, x4, x3, x1, x0);
}

void Task68::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: Gamma395
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x0, x3, x2, x1);
  {
    if (x2 == x3 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x1) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x1) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix1+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x4 && x0 == x3) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                  += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x3) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix3+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (2.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix2+x2.size()*(ix1))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    if (x2 == x3 && x0 == x4) {
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
            for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
              for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                  += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1))];
              }
            }
          }
        }
      }
    }
  }
  {
    if (x0 == x4) {
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x3, x2, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
          for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix4+x0.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix3+x3.size()*(ix2+x2.size()*(ix1))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x4, x0, x1);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix3+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix1))))];
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
      std::unique_ptr<double[]> i0data = in(1)->get_block(ci0, x5, x1, x0, x3);
      for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
        for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
          for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
            for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
              for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0+ci0.size()*(ix5+x5.size()*(ix4+x4.size()*(ix0+x0.size()*(ix3+x3.size()*(ix4+x2.size()*(ix1))))))]
                    += (-1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix1+x1.size()*(ix0+x0.size()*(ix3))))];
                }
              }
            }
          }
        }
      }
    }
  }
  {
    std::unique_ptr<double[]> i0data = in(2)->get_block(ci0, x5, x4, x0, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,1>(i0data, odata, ci0.size(), x5.size(), x4.size(), x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x0, x3, x2, x1);
}

void Task70::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c1, c2, x0, x1");
  madness::World::get_default().gop.fence();
}

void Task71::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x5, x1, x4") * (*ta2_)("c1, x5, c2, x4") * 2;
  madness::World::get_default().gop.fence();
}

void Task72::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c2, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c1, x3, c2, x2");
  madness::World::get_default().gop.fence();
}

void Task73::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x3, c2, x2") += (*ta1_)("c1, x3, c2, x2") * e0_ * (-2)
     + (*ta2_)("c1, x3, c2, x2");
  madness::World::get_default().gop.fence();
}

void Task74::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, x1, c1, x0") += (*ta1_)("c2, c1, x0, x1") + (*ta1_)("c1, c2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task75::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x0, x3, x1, x2") * (*ta2_)("c2, c1, x3, x2");
  madness::World::get_default().gop.fence();
}

void Task76::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c1, x3, x2") += (*ta1_)("c1, x3, c3, x2") * (*ta2_)("c2, c3") * (-2);
  madness::World::get_default().gop.fence();
}

void Task77::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("c2, x2") * (*ta2_)("c1, x0, x1, x2");
  madness::World::get_default().gop.fence();
}

void Task78::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x0, x1, x2") += (*ta1_)("x5, x4, x0, x3, x1, x2") * (*ta2_)("x5, x4, c1, x3");
  madness::World::get_default().gop.fence();
}

void Task79::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c2, c1, x0, x1") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("x2, c1, c2, x3");
  madness::World::get_default().gop.fence();
}

void Task80::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x2, c1, c2, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, x2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task81::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x0, x1") += (*ta1_)("c1, x2, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task82::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x3, x4, x1, x0") * (*ta2_)("x3, c1, x5, x4");
  madness::World::get_default().gop.fence();
}

void Task83::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, c1, x5, x4") += (*ta1_)("c1, x5, c2, x4") * (*ta2_)("x3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task84::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x7, x6, x2, x5, x1, x0") * (*ta2_)("x7, x6, c1, x5");
  madness::World::get_default().gop.fence();
}

void Task85::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x4, x2, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3");
  madness::World::get_default().gop.fence();
}

void Task86::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, c1, x3") * e0_ * (-1)
     + (*ta2_)("x5, x4, c1, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task87::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("x5, x4, c2, x3") * (*ta2_)("c1, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task88::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x5, x4, x3") += (*ta1_)("c1, a2, x5, x4") * (*ta2_)("a2, x3");
  madness::World::get_default().gop.fence();
}

void Task89::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x3, x1, x0") * (*ta2_)("c1, x3");
  madness::World::get_default().gop.fence();
}

void Task90::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x3") += (*ta1_)("c1, x3");
  madness::World::get_default().gop.fence();
}

void Task91::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x3") += (*ta1_)("c2, a3, c1, x3") * (*ta2_)("a3, c2") * 2;
  madness::World::get_default().gop.fence();
}

void Task92::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x3") += (*ta1_)("c1, a3, c2, x3") * (*ta2_)("a3, c2") * (-1);
  madness::World::get_default().gop.fence();
}

void Task93::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x5, x3, x2, x4, x1, x0") * (*ta2_)("x3, x5, c1, x4");
  madness::World::get_default().gop.fence();
}

void Task94::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("x3, x5, c1, x4") += (*ta1_)("x5, a2, c1, x4") * (*ta2_)("a2, x3") * (-1);
  madness::World::get_default().gop.fence();
}

void Task95::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, x2, x1, x0") += (*ta1_)("x2, x5, x4, x3, x1, x0") * (*ta2_)("c1, x5, x4, x3") * 0.5;
  madness::World::get_default().gop.fence();
}

void Task96::compute_() {
  madness::World::get_default().gop.fence();
  (*ta0_)("c3, x0, c1, a2") += (*ta1_)("c1, c3, x0, a2");
  madness::World::get_default().gop.fence();
}

void Task97::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("x1, a2") * (*ta2_)("c1, c3, x1, x0");
  madness::World::get_default().gop.fence();
}

void Task98::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  ta1_->init();
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c3, x1, x0") += (*ta1_)("x1, x3, x0, x2") * (*ta2_)("c1, x3, c3, x2") * (-2);
  madness::World::get_default().gop.fence();
}

void Task99::compute_() {
  if (!ta0_->initialized())
    ta0_->fill_local(0.0);
  madness::World::get_default().gop.fence();
  (*ta0_)("c1, c3, x0, a2") += (*ta1_)("c1, a2") * (*ta2_)("c3, x0");
  madness::World::get_default().gop.fence();
}

#endif
