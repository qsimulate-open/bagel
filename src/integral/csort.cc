//
// BAGEL - Parallel electron correlation program.
// Filename: csort.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <algorithm>
#include <src/integral/sortlist.h>

using namespace std;
using namespace bagel;

void CSortList::sort_indices_00(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 1;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 1 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 1;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 1 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 1 * (c3 + c3end * c2);
          const int toffset = 1 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_01(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 3;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 3 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 3;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 3 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 3 * (c3 + c3end * c2);
          const int toffset = 3 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_11(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 9;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 3 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 3;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 9 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 9 * (c3 + c3end * c2);
          const int toffset = 3 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_02(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 6;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 6 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 6;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 6 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 3];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 4];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 5];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 6 * (c3 + c3end * c2);
          const int toffset = 6 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  3,   1, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+  4,   1, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+  5,   1, current_target+toffset+ 5*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_12(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 18;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 6 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 6;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 18 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 9];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 10];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 17];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 18 * (c3 + c3end * c2);
          const int toffset = 6 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  9,   3, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 12,   3, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 15,   3, current_target+toffset+ 5*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_22(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 36;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 6 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 6;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 36 * (c3 + c3end * c2);
          const int toffset = 6 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 8];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 9];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 14];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 15];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 16];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 20];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 21];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 22];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 26];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 27];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 28];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 35];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 6 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 6;
          const int soffset = 36 * (c3 + c3end * c2);
          const int toffset = 6 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   6, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  6,   6, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 12,   6, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 18,   6, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 24,   6, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 30,   6, current_target+toffset+ 5*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_03(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 10;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 10 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 10;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 10 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 3];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 4];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 6];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 7];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 9];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 10 * (c3 + c3end * c2);
          const int toffset = 10 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  3,   1, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+  4,   1, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+  5,   1, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+  6,   1, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+  7,   1, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+  8,   1, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+  9,   1, current_target+toffset+ 9*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_13(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 30;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 10 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 10;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 30 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 9];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 10];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 26];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 27];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 28];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 29];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 30 * (c3 + c3end * c2);
          const int toffset = 10 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  9,   3, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 12,   3, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 15,   3, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 18,   3, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 21,   3, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 24,   3, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 27,   3, current_target+toffset+ 9*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_23(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 60;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 10 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 10;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 60 * (c3 + c3end * c2);
          const int toffset = 6 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 8];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 9];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 14];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 15];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 16];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 20];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 21];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 22];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 26];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 27];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 28];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 38];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 39];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 40];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 50];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 51];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 52];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 56];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 57];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 58];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 59];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 6 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 6;
          const int soffset = 60 * (c3 + c3end * c2);
          const int toffset = 10 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   6, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  6,   6, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 12,   6, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 18,   6, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 24,   6, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 30,   6, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 36,   6, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 42,   6, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 48,   6, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 54,   6, current_target+toffset+ 9*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_33(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 100;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 10 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 10;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 100 * (c3 + c3end * c2);
          const int toffset = 10 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 12];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 13];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 14];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 20];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 21];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 22];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 23];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 24];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 25];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 26];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 27];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 28];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 39];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 40];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 41];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 42];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 43];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 44];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 45];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 46];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 47];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 48];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 49];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 50];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 51];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 52];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 53];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 54];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 55];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 56];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 57];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 58];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 69];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 70];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 71];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 72];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 73];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 74];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 75];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 76];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 77];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 78];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 79];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 80];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 81];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 82];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 83];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 84];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 85];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 86];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 87];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 88];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 99];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 10 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 10;
          const int soffset = 100 * (c3 + c3end * c2);
          const int toffset = 10 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  10, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 10,  10, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 20,  10, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 30,  10, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 40,  10, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 50,  10, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 60,  10, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 70,  10, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 80,  10, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 90,  10, current_target+toffset+ 9*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_04(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 15;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 15 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 15;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 15 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 3];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 4];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 6];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 7];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 10];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 12];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 13];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 14];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 15 * (c3 + c3end * c2);
          const int toffset = 15 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  3,   1, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+  4,   1, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+  5,   1, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+  6,   1, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+  7,   1, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+  8,   1, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+  9,   1, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 10,   1, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 11,   1, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 12,   1, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 13,   1, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 14,   1, current_target+toffset+14*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_14(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 45;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 15 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 15;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 45 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 9];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 10];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 26];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 27];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 28];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 32];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 33];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 34];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 38];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 39];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 40];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 44];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 45 * (c3 + c3end * c2);
          const int toffset = 15 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  9,   3, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 12,   3, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 15,   3, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 18,   3, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 21,   3, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 24,   3, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 27,   3, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 30,   3, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 33,   3, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 36,   3, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 39,   3, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 42,   3, current_target+toffset+14*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_24(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 90;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 15 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 15;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 90 * (c3 + c3end * c2);
          const int toffset = 6 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 8];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 9];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 14];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 15];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 16];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 20];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 21];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 22];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 26];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 27];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 28];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 38];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 39];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 40];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 50];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 51];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 52];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 56];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 57];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 58];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 65];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 66];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 67];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 68];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 69];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 70];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 71];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 72];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 73];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 74];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 75];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 76];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 77];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 78];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 79];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 80];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 81];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 82];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 89];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 6 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 6;
          const int soffset = 90 * (c3 + c3end * c2);
          const int toffset = 15 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   6, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  6,   6, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 12,   6, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 18,   6, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 24,   6, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 30,   6, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 36,   6, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 42,   6, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 48,   6, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 54,   6, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 60,   6, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 66,   6, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 72,   6, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 78,   6, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 84,   6, current_target+toffset+14*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_34(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 150;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 15 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 15;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 150 * (c3 + c3end * c2);
          const int toffset = 10 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 12];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 13];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 14];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 20];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 21];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 22];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 23];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 24];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 25];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 26];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 27];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 28];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 39];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 40];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 41];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 42];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 43];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 44];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 45];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 46];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 47];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 48];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 49];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 50];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 51];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 52];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 53];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 54];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 55];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 56];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 57];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 58];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 69];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 70];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 71];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 72];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 73];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 74];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 75];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 76];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 77];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 78];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 79];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 80];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 81];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 82];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 83];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 84];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 85];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 86];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 87];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 88];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 99];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 100];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 101];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 102];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 103];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 104];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 105];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 106];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 107];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 108];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 109];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 110];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 111];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 112];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 113];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 114];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 115];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 116];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 117];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 118];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 129];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 130];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 131];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 132];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 133];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 134];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 135];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 136];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 137];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 138];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 139];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 140];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 141];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 142];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 143];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 144];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 145];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 146];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 147];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 148];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 149];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 10 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 10;
          const int soffset = 150 * (c3 + c3end * c2);
          const int toffset = 15 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  10, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 10,  10, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 20,  10, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 30,  10, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 40,  10, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 50,  10, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 60,  10, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 70,  10, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 80,  10, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 90,  10, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+100,  10, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+110,  10, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+120,  10, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+130,  10, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+140,  10, current_target+toffset+14*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_44(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 225;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 15 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 15;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 225 * (c3 + c3end * c2);
          const int toffset = 15 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 20];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 21];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 22];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 23];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 24];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 25];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 26];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 27];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 39];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 40];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 41];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 42];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 43];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 44];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 45];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 46];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 47];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 48];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 49];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 50];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 51];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 52];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 53];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 54];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 55];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 56];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 57];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 58];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 69];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 70];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 71];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 72];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 73];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 74];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 75];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 76];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 77];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 78];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 79];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 80];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 81];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 82];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 83];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 84];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 85];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 86];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 87];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 88];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 99];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 100];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 101];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 102];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 103];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 104];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 105];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 106];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 107];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 108];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 109];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 110];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 111];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 112];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 113];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 114];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 115];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 116];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 117];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 118];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 129];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 130];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 131];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 132];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 133];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 134];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 135];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 136];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 137];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 138];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 139];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 140];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 141];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 142];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 143];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 144];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 145];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 146];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 147];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 148];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 155];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 156];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 157];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 158];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 159];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 160];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 161];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 162];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 163];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 164];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 165];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 166];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 167];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 168];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 169];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 170];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 171];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 172];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 173];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 174];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 175];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 176];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 177];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 178];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 179];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 180];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 181];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 182];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 183];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 184];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 185];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 186];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 187];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 188];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 189];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 190];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 191];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 192];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 193];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 194];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 195];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 196];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 197];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 198];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 199];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 200];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 201];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 202];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 203];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 204];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 205];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 206];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 207];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 208];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 219];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 220];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 221];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 222];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 223];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 224];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 15 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 15;
          const int soffset = 225 * (c3 + c3end * c2);
          const int toffset = 15 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  15, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 15,  15, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 30,  15, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 45,  15, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 60,  15, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 75,  15, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 90,  15, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+105,  15, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+120,  15, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+135,  15, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+150,  15, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+165,  15, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+180,  15, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+195,  15, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+210,  15, current_target+toffset+14*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_05(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 21;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 21 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 3];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 4];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 6];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 7];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 10];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 12];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 13];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 15];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 16];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 18];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 20];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 21 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  3,   1, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+  4,   1, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+  5,   1, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+  6,   1, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+  7,   1, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+  8,   1, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+  9,   1, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 10,   1, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 11,   1, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 12,   1, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 13,   1, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 14,   1, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 15,   1, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 16,   1, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+ 17,   1, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+ 18,   1, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+ 19,   1, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+ 20,   1, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_15(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 63;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 63 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 9];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 10];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 26];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 27];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 28];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 32];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 33];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 34];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 38];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 39];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 40];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 44];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 45];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 46];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 50];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 51];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 52];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 56];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 57];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 58];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 62];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 63 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  9,   3, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 12,   3, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 15,   3, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 18,   3, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 21,   3, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 24,   3, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 27,   3, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 30,   3, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 33,   3, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 36,   3, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 39,   3, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 42,   3, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 45,   3, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 48,   3, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+ 51,   3, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+ 54,   3, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+ 57,   3, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+ 60,   3, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_25(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 126;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 126 * (c3 + c3end * c2);
          const int toffset = 6 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 8];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 9];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 14];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 15];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 16];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 20];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 21];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 22];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 26];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 27];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 28];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 38];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 39];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 40];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 50];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 51];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 52];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 56];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 57];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 58];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 65];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 66];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 67];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 68];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 69];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 70];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 71];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 72];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 73];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 74];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 75];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 76];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 77];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 78];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 79];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 80];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 81];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 82];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 95];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 96];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 97];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 98];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 99];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 100];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 101];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 102];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 103];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 104];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 105];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 106];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 107];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 108];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 109];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 110];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 111];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 112];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 113];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 114];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 115];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 116];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 117];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 118];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 125];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 6 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 6;
          const int soffset = 126 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   6, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  6,   6, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 12,   6, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 18,   6, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 24,   6, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 30,   6, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 36,   6, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 42,   6, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 48,   6, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 54,   6, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 60,   6, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 66,   6, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 72,   6, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 78,   6, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 84,   6, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 90,   6, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 96,   6, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+102,   6, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+108,   6, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+114,   6, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+120,   6, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_35(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 210;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 210 * (c3 + c3end * c2);
          const int toffset = 10 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 12];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 13];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 14];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 20];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 21];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 22];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 23];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 24];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 25];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 26];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 27];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 28];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 39];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 40];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 41];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 42];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 43];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 44];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 45];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 46];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 47];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 48];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 49];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 50];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 51];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 52];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 53];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 54];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 55];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 56];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 57];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 58];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 69];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 70];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 71];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 72];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 73];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 74];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 75];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 76];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 77];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 78];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 79];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 80];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 81];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 82];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 83];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 84];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 85];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 86];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 87];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 88];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 99];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 100];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 101];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 102];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 103];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 104];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 105];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 106];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 107];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 108];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 109];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 110];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 111];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 112];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 113];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 114];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 115];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 116];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 117];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 118];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 129];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 130];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 131];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 132];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 133];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 134];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 135];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 136];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 137];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 138];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 139];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 140];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 141];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 142];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 143];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 144];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 145];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 146];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 147];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 148];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 155];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 156];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 157];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 158];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 159];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 160];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 161];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 162];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 163];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 164];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 165];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 166];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 167];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 168];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 169];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 170];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 171];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 172];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 173];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 174];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 175];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 176];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 177];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 178];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 179];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 180];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 181];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 182];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 183];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 184];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 185];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 186];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 187];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 188];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 189];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 190];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 191];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 192];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 193];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 194];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 195];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 196];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 197];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 198];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 199];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 200];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 201];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 202];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 203];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 204];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 205];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 206];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 207];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 208];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 209];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 10 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 10;
          const int soffset = 210 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  10, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 10,  10, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 20,  10, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 30,  10, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 40,  10, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 50,  10, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 60,  10, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 70,  10, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 80,  10, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 90,  10, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+100,  10, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+110,  10, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+120,  10, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+130,  10, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+140,  10, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+150,  10, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+160,  10, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+170,  10, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+180,  10, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+190,  10, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+200,  10, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_45(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 315;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 315 * (c3 + c3end * c2);
          const int toffset = 15 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 20];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 21];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 22];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 23];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 24];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 25];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 26];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 27];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 39];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 40];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 41];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 42];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 43];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 44];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 45];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 46];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 47];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 48];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 49];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 50];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 51];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 52];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 53];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 54];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 55];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 56];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 57];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 58];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 69];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 70];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 71];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 72];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 73];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 74];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 75];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 76];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 77];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 78];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 79];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 80];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 81];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 82];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 83];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 84];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 85];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 86];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 87];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 88];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 99];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 100];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 101];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 102];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 103];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 104];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 105];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 106];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 107];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 108];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 109];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 110];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 111];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 112];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 113];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 114];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 115];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 116];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 117];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 118];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 129];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 130];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 131];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 132];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 133];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 134];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 135];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 136];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 137];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 138];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 139];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 140];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 141];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 142];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 143];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 144];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 145];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 146];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 147];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 148];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 155];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 156];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 157];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 158];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 159];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 160];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 161];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 162];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 163];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 164];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 165];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 166];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 167];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 168];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 169];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 170];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 171];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 172];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 173];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 174];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 175];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 176];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 177];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 178];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 179];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 180];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 181];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 182];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 183];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 184];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 185];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 186];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 187];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 188];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 189];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 190];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 191];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 192];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 193];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 194];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 195];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 196];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 197];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 198];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 199];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 200];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 201];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 202];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 203];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 204];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 205];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 206];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 207];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 208];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 219];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 220];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 221];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 222];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 223];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 224];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 225];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 226];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 227];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 228];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 229];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 230];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 231];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 232];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 233];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 234];
          current_target[toffset + 10 * cont2csize + 15] = current_source[soffset + 235];
          current_target[toffset + 11 * cont2csize + 15] = current_source[soffset + 236];
          current_target[toffset + 12 * cont2csize + 15] = current_source[soffset + 237];
          current_target[toffset + 13 * cont2csize + 15] = current_source[soffset + 238];
          current_target[toffset + 14 * cont2csize + 15] = current_source[soffset + 239];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 240];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 241];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 242];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 243];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 244];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 245];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 246];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 247];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 248];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 249];
          current_target[toffset + 10 * cont2csize + 16] = current_source[soffset + 250];
          current_target[toffset + 11 * cont2csize + 16] = current_source[soffset + 251];
          current_target[toffset + 12 * cont2csize + 16] = current_source[soffset + 252];
          current_target[toffset + 13 * cont2csize + 16] = current_source[soffset + 253];
          current_target[toffset + 14 * cont2csize + 16] = current_source[soffset + 254];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 255];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 256];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 257];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 258];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 259];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 260];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 261];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 262];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 263];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 264];
          current_target[toffset + 10 * cont2csize + 17] = current_source[soffset + 265];
          current_target[toffset + 11 * cont2csize + 17] = current_source[soffset + 266];
          current_target[toffset + 12 * cont2csize + 17] = current_source[soffset + 267];
          current_target[toffset + 13 * cont2csize + 17] = current_source[soffset + 268];
          current_target[toffset + 14 * cont2csize + 17] = current_source[soffset + 269];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 270];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 271];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 272];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 273];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 274];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 275];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 276];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 277];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 278];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 279];
          current_target[toffset + 10 * cont2csize + 18] = current_source[soffset + 280];
          current_target[toffset + 11 * cont2csize + 18] = current_source[soffset + 281];
          current_target[toffset + 12 * cont2csize + 18] = current_source[soffset + 282];
          current_target[toffset + 13 * cont2csize + 18] = current_source[soffset + 283];
          current_target[toffset + 14 * cont2csize + 18] = current_source[soffset + 284];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 285];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 286];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 287];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 288];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 289];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 290];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 291];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 292];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 293];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 294];
          current_target[toffset + 10 * cont2csize + 19] = current_source[soffset + 295];
          current_target[toffset + 11 * cont2csize + 19] = current_source[soffset + 296];
          current_target[toffset + 12 * cont2csize + 19] = current_source[soffset + 297];
          current_target[toffset + 13 * cont2csize + 19] = current_source[soffset + 298];
          current_target[toffset + 14 * cont2csize + 19] = current_source[soffset + 299];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 300];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 301];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 302];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 303];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 304];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 305];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 306];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 307];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 308];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 309];
          current_target[toffset + 10 * cont2csize + 20] = current_source[soffset + 310];
          current_target[toffset + 11 * cont2csize + 20] = current_source[soffset + 311];
          current_target[toffset + 12 * cont2csize + 20] = current_source[soffset + 312];
          current_target[toffset + 13 * cont2csize + 20] = current_source[soffset + 313];
          current_target[toffset + 14 * cont2csize + 20] = current_source[soffset + 314];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 15 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 15;
          const int soffset = 315 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  15, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 15,  15, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 30,  15, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 45,  15, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 60,  15, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 75,  15, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 90,  15, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+105,  15, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+120,  15, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+135,  15, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+150,  15, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+165,  15, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+180,  15, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+195,  15, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+210,  15, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+225,  15, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+240,  15, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+255,  15, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+270,  15, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+285,  15, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+300,  15, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_55(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 441;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 21 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 21;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 441 * (c3 + c3end * c2);
          const int toffset = 21 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 15 * cont2csize + 0] = current_source[soffset + 15];
          current_target[toffset + 16 * cont2csize + 0] = current_source[soffset + 16];
          current_target[toffset + 17 * cont2csize + 0] = current_source[soffset + 17];
          current_target[toffset + 18 * cont2csize + 0] = current_source[soffset + 18];
          current_target[toffset + 19 * cont2csize + 0] = current_source[soffset + 19];
          current_target[toffset + 20 * cont2csize + 0] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 23];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 24];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 25];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 26];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 27];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 30];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 31];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 32];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 33];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 34];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 35];
          current_target[toffset + 15 * cont2csize + 1] = current_source[soffset + 36];
          current_target[toffset + 16 * cont2csize + 1] = current_source[soffset + 37];
          current_target[toffset + 17 * cont2csize + 1] = current_source[soffset + 38];
          current_target[toffset + 18 * cont2csize + 1] = current_source[soffset + 39];
          current_target[toffset + 19 * cont2csize + 1] = current_source[soffset + 40];
          current_target[toffset + 20 * cont2csize + 1] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 47];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 48];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 49];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 50];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 51];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 52];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 53];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 54];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 55];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 56];
          current_target[toffset + 15 * cont2csize + 2] = current_source[soffset + 57];
          current_target[toffset + 16 * cont2csize + 2] = current_source[soffset + 58];
          current_target[toffset + 17 * cont2csize + 2] = current_source[soffset + 59];
          current_target[toffset + 18 * cont2csize + 2] = current_source[soffset + 60];
          current_target[toffset + 19 * cont2csize + 2] = current_source[soffset + 61];
          current_target[toffset + 20 * cont2csize + 2] = current_source[soffset + 62];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 63];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 64];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 65];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 66];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 67];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 68];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 69];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 70];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 71];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 72];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 73];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 74];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 75];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 76];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 77];
          current_target[toffset + 15 * cont2csize + 3] = current_source[soffset + 78];
          current_target[toffset + 16 * cont2csize + 3] = current_source[soffset + 79];
          current_target[toffset + 17 * cont2csize + 3] = current_source[soffset + 80];
          current_target[toffset + 18 * cont2csize + 3] = current_source[soffset + 81];
          current_target[toffset + 19 * cont2csize + 3] = current_source[soffset + 82];
          current_target[toffset + 20 * cont2csize + 3] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 89];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 90];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 91];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 92];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 93];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 94];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 95];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 96];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 97];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 98];
          current_target[toffset + 15 * cont2csize + 4] = current_source[soffset + 99];
          current_target[toffset + 16 * cont2csize + 4] = current_source[soffset + 100];
          current_target[toffset + 17 * cont2csize + 4] = current_source[soffset + 101];
          current_target[toffset + 18 * cont2csize + 4] = current_source[soffset + 102];
          current_target[toffset + 19 * cont2csize + 4] = current_source[soffset + 103];
          current_target[toffset + 20 * cont2csize + 4] = current_source[soffset + 104];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 105];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 106];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 107];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 108];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 109];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 110];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 111];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 112];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 113];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 114];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 115];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 116];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 117];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 118];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 119];
          current_target[toffset + 15 * cont2csize + 5] = current_source[soffset + 120];
          current_target[toffset + 16 * cont2csize + 5] = current_source[soffset + 121];
          current_target[toffset + 17 * cont2csize + 5] = current_source[soffset + 122];
          current_target[toffset + 18 * cont2csize + 5] = current_source[soffset + 123];
          current_target[toffset + 19 * cont2csize + 5] = current_source[soffset + 124];
          current_target[toffset + 20 * cont2csize + 5] = current_source[soffset + 125];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 126];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 127];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 128];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 129];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 130];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 131];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 132];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 133];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 134];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 135];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 136];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 137];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 138];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 139];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 140];
          current_target[toffset + 15 * cont2csize + 6] = current_source[soffset + 141];
          current_target[toffset + 16 * cont2csize + 6] = current_source[soffset + 142];
          current_target[toffset + 17 * cont2csize + 6] = current_source[soffset + 143];
          current_target[toffset + 18 * cont2csize + 6] = current_source[soffset + 144];
          current_target[toffset + 19 * cont2csize + 6] = current_source[soffset + 145];
          current_target[toffset + 20 * cont2csize + 6] = current_source[soffset + 146];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 147];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 148];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 149];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 150];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 151];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 152];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 153];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 154];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 155];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 156];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 157];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 158];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 159];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 160];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 161];
          current_target[toffset + 15 * cont2csize + 7] = current_source[soffset + 162];
          current_target[toffset + 16 * cont2csize + 7] = current_source[soffset + 163];
          current_target[toffset + 17 * cont2csize + 7] = current_source[soffset + 164];
          current_target[toffset + 18 * cont2csize + 7] = current_source[soffset + 165];
          current_target[toffset + 19 * cont2csize + 7] = current_source[soffset + 166];
          current_target[toffset + 20 * cont2csize + 7] = current_source[soffset + 167];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 168];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 169];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 170];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 171];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 172];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 173];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 174];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 175];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 176];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 177];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 178];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 179];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 180];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 181];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 182];
          current_target[toffset + 15 * cont2csize + 8] = current_source[soffset + 183];
          current_target[toffset + 16 * cont2csize + 8] = current_source[soffset + 184];
          current_target[toffset + 17 * cont2csize + 8] = current_source[soffset + 185];
          current_target[toffset + 18 * cont2csize + 8] = current_source[soffset + 186];
          current_target[toffset + 19 * cont2csize + 8] = current_source[soffset + 187];
          current_target[toffset + 20 * cont2csize + 8] = current_source[soffset + 188];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 189];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 190];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 191];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 192];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 193];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 194];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 195];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 196];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 197];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 198];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 199];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 200];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 201];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 202];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 203];
          current_target[toffset + 15 * cont2csize + 9] = current_source[soffset + 204];
          current_target[toffset + 16 * cont2csize + 9] = current_source[soffset + 205];
          current_target[toffset + 17 * cont2csize + 9] = current_source[soffset + 206];
          current_target[toffset + 18 * cont2csize + 9] = current_source[soffset + 207];
          current_target[toffset + 19 * cont2csize + 9] = current_source[soffset + 208];
          current_target[toffset + 20 * cont2csize + 9] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 219];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 220];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 221];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 222];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 223];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 224];
          current_target[toffset + 15 * cont2csize + 10] = current_source[soffset + 225];
          current_target[toffset + 16 * cont2csize + 10] = current_source[soffset + 226];
          current_target[toffset + 17 * cont2csize + 10] = current_source[soffset + 227];
          current_target[toffset + 18 * cont2csize + 10] = current_source[soffset + 228];
          current_target[toffset + 19 * cont2csize + 10] = current_source[soffset + 229];
          current_target[toffset + 20 * cont2csize + 10] = current_source[soffset + 230];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 231];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 232];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 233];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 234];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 235];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 236];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 237];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 238];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 239];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 240];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 241];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 242];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 243];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 244];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 245];
          current_target[toffset + 15 * cont2csize + 11] = current_source[soffset + 246];
          current_target[toffset + 16 * cont2csize + 11] = current_source[soffset + 247];
          current_target[toffset + 17 * cont2csize + 11] = current_source[soffset + 248];
          current_target[toffset + 18 * cont2csize + 11] = current_source[soffset + 249];
          current_target[toffset + 19 * cont2csize + 11] = current_source[soffset + 250];
          current_target[toffset + 20 * cont2csize + 11] = current_source[soffset + 251];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 252];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 253];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 254];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 255];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 256];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 257];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 258];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 259];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 260];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 261];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 262];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 263];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 264];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 265];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 266];
          current_target[toffset + 15 * cont2csize + 12] = current_source[soffset + 267];
          current_target[toffset + 16 * cont2csize + 12] = current_source[soffset + 268];
          current_target[toffset + 17 * cont2csize + 12] = current_source[soffset + 269];
          current_target[toffset + 18 * cont2csize + 12] = current_source[soffset + 270];
          current_target[toffset + 19 * cont2csize + 12] = current_source[soffset + 271];
          current_target[toffset + 20 * cont2csize + 12] = current_source[soffset + 272];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 273];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 274];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 275];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 276];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 277];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 278];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 279];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 280];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 281];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 282];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 283];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 284];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 285];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 286];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 287];
          current_target[toffset + 15 * cont2csize + 13] = current_source[soffset + 288];
          current_target[toffset + 16 * cont2csize + 13] = current_source[soffset + 289];
          current_target[toffset + 17 * cont2csize + 13] = current_source[soffset + 290];
          current_target[toffset + 18 * cont2csize + 13] = current_source[soffset + 291];
          current_target[toffset + 19 * cont2csize + 13] = current_source[soffset + 292];
          current_target[toffset + 20 * cont2csize + 13] = current_source[soffset + 293];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 294];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 295];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 296];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 297];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 298];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 299];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 300];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 301];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 302];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 303];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 304];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 305];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 306];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 307];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 308];
          current_target[toffset + 15 * cont2csize + 14] = current_source[soffset + 309];
          current_target[toffset + 16 * cont2csize + 14] = current_source[soffset + 310];
          current_target[toffset + 17 * cont2csize + 14] = current_source[soffset + 311];
          current_target[toffset + 18 * cont2csize + 14] = current_source[soffset + 312];
          current_target[toffset + 19 * cont2csize + 14] = current_source[soffset + 313];
          current_target[toffset + 20 * cont2csize + 14] = current_source[soffset + 314];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 315];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 316];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 317];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 318];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 319];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 320];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 321];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 322];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 323];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 324];
          current_target[toffset + 10 * cont2csize + 15] = current_source[soffset + 325];
          current_target[toffset + 11 * cont2csize + 15] = current_source[soffset + 326];
          current_target[toffset + 12 * cont2csize + 15] = current_source[soffset + 327];
          current_target[toffset + 13 * cont2csize + 15] = current_source[soffset + 328];
          current_target[toffset + 14 * cont2csize + 15] = current_source[soffset + 329];
          current_target[toffset + 15 * cont2csize + 15] = current_source[soffset + 330];
          current_target[toffset + 16 * cont2csize + 15] = current_source[soffset + 331];
          current_target[toffset + 17 * cont2csize + 15] = current_source[soffset + 332];
          current_target[toffset + 18 * cont2csize + 15] = current_source[soffset + 333];
          current_target[toffset + 19 * cont2csize + 15] = current_source[soffset + 334];
          current_target[toffset + 20 * cont2csize + 15] = current_source[soffset + 335];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 336];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 337];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 338];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 339];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 340];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 341];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 342];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 343];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 344];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 345];
          current_target[toffset + 10 * cont2csize + 16] = current_source[soffset + 346];
          current_target[toffset + 11 * cont2csize + 16] = current_source[soffset + 347];
          current_target[toffset + 12 * cont2csize + 16] = current_source[soffset + 348];
          current_target[toffset + 13 * cont2csize + 16] = current_source[soffset + 349];
          current_target[toffset + 14 * cont2csize + 16] = current_source[soffset + 350];
          current_target[toffset + 15 * cont2csize + 16] = current_source[soffset + 351];
          current_target[toffset + 16 * cont2csize + 16] = current_source[soffset + 352];
          current_target[toffset + 17 * cont2csize + 16] = current_source[soffset + 353];
          current_target[toffset + 18 * cont2csize + 16] = current_source[soffset + 354];
          current_target[toffset + 19 * cont2csize + 16] = current_source[soffset + 355];
          current_target[toffset + 20 * cont2csize + 16] = current_source[soffset + 356];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 357];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 358];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 359];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 360];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 361];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 362];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 363];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 364];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 365];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 366];
          current_target[toffset + 10 * cont2csize + 17] = current_source[soffset + 367];
          current_target[toffset + 11 * cont2csize + 17] = current_source[soffset + 368];
          current_target[toffset + 12 * cont2csize + 17] = current_source[soffset + 369];
          current_target[toffset + 13 * cont2csize + 17] = current_source[soffset + 370];
          current_target[toffset + 14 * cont2csize + 17] = current_source[soffset + 371];
          current_target[toffset + 15 * cont2csize + 17] = current_source[soffset + 372];
          current_target[toffset + 16 * cont2csize + 17] = current_source[soffset + 373];
          current_target[toffset + 17 * cont2csize + 17] = current_source[soffset + 374];
          current_target[toffset + 18 * cont2csize + 17] = current_source[soffset + 375];
          current_target[toffset + 19 * cont2csize + 17] = current_source[soffset + 376];
          current_target[toffset + 20 * cont2csize + 17] = current_source[soffset + 377];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 378];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 379];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 380];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 381];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 382];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 383];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 384];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 385];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 386];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 387];
          current_target[toffset + 10 * cont2csize + 18] = current_source[soffset + 388];
          current_target[toffset + 11 * cont2csize + 18] = current_source[soffset + 389];
          current_target[toffset + 12 * cont2csize + 18] = current_source[soffset + 390];
          current_target[toffset + 13 * cont2csize + 18] = current_source[soffset + 391];
          current_target[toffset + 14 * cont2csize + 18] = current_source[soffset + 392];
          current_target[toffset + 15 * cont2csize + 18] = current_source[soffset + 393];
          current_target[toffset + 16 * cont2csize + 18] = current_source[soffset + 394];
          current_target[toffset + 17 * cont2csize + 18] = current_source[soffset + 395];
          current_target[toffset + 18 * cont2csize + 18] = current_source[soffset + 396];
          current_target[toffset + 19 * cont2csize + 18] = current_source[soffset + 397];
          current_target[toffset + 20 * cont2csize + 18] = current_source[soffset + 398];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 399];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 400];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 401];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 402];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 403];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 404];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 405];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 406];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 407];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 408];
          current_target[toffset + 10 * cont2csize + 19] = current_source[soffset + 409];
          current_target[toffset + 11 * cont2csize + 19] = current_source[soffset + 410];
          current_target[toffset + 12 * cont2csize + 19] = current_source[soffset + 411];
          current_target[toffset + 13 * cont2csize + 19] = current_source[soffset + 412];
          current_target[toffset + 14 * cont2csize + 19] = current_source[soffset + 413];
          current_target[toffset + 15 * cont2csize + 19] = current_source[soffset + 414];
          current_target[toffset + 16 * cont2csize + 19] = current_source[soffset + 415];
          current_target[toffset + 17 * cont2csize + 19] = current_source[soffset + 416];
          current_target[toffset + 18 * cont2csize + 19] = current_source[soffset + 417];
          current_target[toffset + 19 * cont2csize + 19] = current_source[soffset + 418];
          current_target[toffset + 20 * cont2csize + 19] = current_source[soffset + 419];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 420];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 421];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 422];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 423];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 424];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 425];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 426];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 427];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 428];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 429];
          current_target[toffset + 10 * cont2csize + 20] = current_source[soffset + 430];
          current_target[toffset + 11 * cont2csize + 20] = current_source[soffset + 431];
          current_target[toffset + 12 * cont2csize + 20] = current_source[soffset + 432];
          current_target[toffset + 13 * cont2csize + 20] = current_source[soffset + 433];
          current_target[toffset + 14 * cont2csize + 20] = current_source[soffset + 434];
          current_target[toffset + 15 * cont2csize + 20] = current_source[soffset + 435];
          current_target[toffset + 16 * cont2csize + 20] = current_source[soffset + 436];
          current_target[toffset + 17 * cont2csize + 20] = current_source[soffset + 437];
          current_target[toffset + 18 * cont2csize + 20] = current_source[soffset + 438];
          current_target[toffset + 19 * cont2csize + 20] = current_source[soffset + 439];
          current_target[toffset + 20 * cont2csize + 20] = current_source[soffset + 440];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 21 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 21;
          const int soffset = 441 * (c3 + c3end * c2);
          const int toffset = 21 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  21, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 21,  21, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 42,  21, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 63,  21, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 84,  21, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+105,  21, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+126,  21, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+147,  21, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+168,  21, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+189,  21, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+210,  21, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+231,  21, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+252,  21, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+273,  21, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+294,  21, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+315,  21, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+336,  21, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+357,  21, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+378,  21, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+399,  21, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+420,  21, current_target+toffset+20*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_06(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 28;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 28 * (c3 + c3end * c2);
          const int toffset = 1 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 1];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 3];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 4];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 6];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 7];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 10];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 12];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 13];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 15];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 16];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 18];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 21];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 22];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 24];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 25];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 26];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 27];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 1 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 1;
          const int soffset = 28 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   1, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  1,   1, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  2,   1, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  3,   1, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+  4,   1, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+  5,   1, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+  6,   1, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+  7,   1, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+  8,   1, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+  9,   1, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 10,   1, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 11,   1, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 12,   1, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 13,   1, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 14,   1, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 15,   1, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 16,   1, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+ 17,   1, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+ 18,   1, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+ 19,   1, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+ 20,   1, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+ 21,   1, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+ 22,   1, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+ 23,   1, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+ 24,   1, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+ 25,   1, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+ 26,   1, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+ 27,   1, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_16(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 84;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 84 * (c3 + c3end * c2);
          const int toffset = 3 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 3];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 4];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 8];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 9];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 10];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 26];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 27];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 28];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 32];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 33];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 34];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 38];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 39];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 40];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 44];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 45];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 46];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 50];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 51];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 52];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 56];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 57];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 58];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 62];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 63];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 64];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 65];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 66];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 67];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 68];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 69];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 70];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 71];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 72];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 73];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 74];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 75];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 76];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 77];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 78];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 79];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 80];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 81];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 82];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 83];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 3 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 3;
          const int soffset = 84 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   3, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  3,   3, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+  6,   3, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+  9,   3, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 12,   3, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 15,   3, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 18,   3, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 21,   3, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 24,   3, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 27,   3, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 30,   3, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 33,   3, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 36,   3, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 39,   3, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 42,   3, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 45,   3, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 48,   3, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+ 51,   3, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+ 54,   3, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+ 57,   3, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+ 60,   3, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+ 63,   3, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+ 66,   3, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+ 69,   3, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+ 72,   3, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+ 75,   3, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+ 78,   3, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+ 81,   3, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_26(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 168;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 168 * (c3 + c3end * c2);
          const int toffset = 6 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 6];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 7];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 8];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 9];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 12];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 13];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 14];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 15];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 16];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 17];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 18];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 19];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 20];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 21];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 22];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 23];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 24];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 25];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 26];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 27];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 28];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 35];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 36];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 37];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 38];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 39];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 40];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 47];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 48];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 49];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 50];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 51];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 52];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 53];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 54];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 55];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 56];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 57];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 58];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 65];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 66];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 67];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 68];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 69];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 70];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 71];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 72];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 73];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 74];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 75];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 76];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 77];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 78];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 79];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 80];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 81];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 82];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 95];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 96];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 97];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 98];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 99];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 100];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 101];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 102];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 103];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 104];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 105];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 106];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 107];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 108];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 109];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 110];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 111];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 112];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 113];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 114];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 115];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 116];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 117];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 118];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 125];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 126];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 127];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 128];
          current_target[toffset + 3 * cont2csize + 21] = current_source[soffset + 129];
          current_target[toffset + 4 * cont2csize + 21] = current_source[soffset + 130];
          current_target[toffset + 5 * cont2csize + 21] = current_source[soffset + 131];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 132];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 133];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 134];
          current_target[toffset + 3 * cont2csize + 22] = current_source[soffset + 135];
          current_target[toffset + 4 * cont2csize + 22] = current_source[soffset + 136];
          current_target[toffset + 5 * cont2csize + 22] = current_source[soffset + 137];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 138];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 139];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 140];
          current_target[toffset + 3 * cont2csize + 23] = current_source[soffset + 141];
          current_target[toffset + 4 * cont2csize + 23] = current_source[soffset + 142];
          current_target[toffset + 5 * cont2csize + 23] = current_source[soffset + 143];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 144];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 145];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 146];
          current_target[toffset + 3 * cont2csize + 24] = current_source[soffset + 147];
          current_target[toffset + 4 * cont2csize + 24] = current_source[soffset + 148];
          current_target[toffset + 5 * cont2csize + 24] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 25] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 25] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 25] = current_source[soffset + 155];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 156];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 157];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 158];
          current_target[toffset + 3 * cont2csize + 26] = current_source[soffset + 159];
          current_target[toffset + 4 * cont2csize + 26] = current_source[soffset + 160];
          current_target[toffset + 5 * cont2csize + 26] = current_source[soffset + 161];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 162];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 163];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 164];
          current_target[toffset + 3 * cont2csize + 27] = current_source[soffset + 165];
          current_target[toffset + 4 * cont2csize + 27] = current_source[soffset + 166];
          current_target[toffset + 5 * cont2csize + 27] = current_source[soffset + 167];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 6 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 6;
          const int soffset = 168 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,   6, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+  6,   6, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 12,   6, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 18,   6, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 24,   6, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 30,   6, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 36,   6, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 42,   6, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 48,   6, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 54,   6, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+ 60,   6, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+ 66,   6, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+ 72,   6, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+ 78,   6, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+ 84,   6, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+ 90,   6, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+ 96,   6, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+102,   6, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+108,   6, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+114,   6, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+120,   6, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+126,   6, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+132,   6, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+138,   6, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+144,   6, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+150,   6, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+156,   6, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+162,   6, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_36(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 280;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 280 * (c3 + c3end * c2);
          const int toffset = 10 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 10];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 11];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 12];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 13];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 14];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 20];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 21];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 22];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 23];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 24];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 25];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 26];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 27];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 28];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 39];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 40];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 41];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 42];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 43];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 44];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 45];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 46];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 47];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 48];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 49];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 50];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 51];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 52];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 53];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 54];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 55];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 56];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 57];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 58];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 69];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 70];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 71];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 72];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 73];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 74];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 75];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 76];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 77];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 78];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 79];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 80];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 81];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 82];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 83];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 84];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 85];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 86];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 87];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 88];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 99];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 100];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 101];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 102];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 103];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 104];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 105];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 106];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 107];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 108];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 109];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 110];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 111];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 112];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 113];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 114];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 115];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 116];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 117];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 118];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 129];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 130];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 131];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 132];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 133];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 134];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 135];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 136];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 137];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 138];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 139];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 140];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 141];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 142];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 143];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 144];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 145];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 146];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 147];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 148];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 155];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 156];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 157];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 158];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 159];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 160];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 161];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 162];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 163];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 164];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 165];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 166];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 167];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 168];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 169];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 170];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 171];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 172];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 173];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 174];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 175];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 176];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 177];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 178];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 179];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 180];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 181];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 182];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 183];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 184];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 185];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 186];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 187];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 188];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 189];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 190];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 191];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 192];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 193];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 194];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 195];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 196];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 197];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 198];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 199];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 200];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 201];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 202];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 203];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 204];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 205];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 206];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 207];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 208];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 21] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 21] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 21] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 21] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 21] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 21] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 21] = current_source[soffset + 219];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 220];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 221];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 222];
          current_target[toffset + 3 * cont2csize + 22] = current_source[soffset + 223];
          current_target[toffset + 4 * cont2csize + 22] = current_source[soffset + 224];
          current_target[toffset + 5 * cont2csize + 22] = current_source[soffset + 225];
          current_target[toffset + 6 * cont2csize + 22] = current_source[soffset + 226];
          current_target[toffset + 7 * cont2csize + 22] = current_source[soffset + 227];
          current_target[toffset + 8 * cont2csize + 22] = current_source[soffset + 228];
          current_target[toffset + 9 * cont2csize + 22] = current_source[soffset + 229];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 230];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 231];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 232];
          current_target[toffset + 3 * cont2csize + 23] = current_source[soffset + 233];
          current_target[toffset + 4 * cont2csize + 23] = current_source[soffset + 234];
          current_target[toffset + 5 * cont2csize + 23] = current_source[soffset + 235];
          current_target[toffset + 6 * cont2csize + 23] = current_source[soffset + 236];
          current_target[toffset + 7 * cont2csize + 23] = current_source[soffset + 237];
          current_target[toffset + 8 * cont2csize + 23] = current_source[soffset + 238];
          current_target[toffset + 9 * cont2csize + 23] = current_source[soffset + 239];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 240];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 241];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 242];
          current_target[toffset + 3 * cont2csize + 24] = current_source[soffset + 243];
          current_target[toffset + 4 * cont2csize + 24] = current_source[soffset + 244];
          current_target[toffset + 5 * cont2csize + 24] = current_source[soffset + 245];
          current_target[toffset + 6 * cont2csize + 24] = current_source[soffset + 246];
          current_target[toffset + 7 * cont2csize + 24] = current_source[soffset + 247];
          current_target[toffset + 8 * cont2csize + 24] = current_source[soffset + 248];
          current_target[toffset + 9 * cont2csize + 24] = current_source[soffset + 249];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 250];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 251];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 252];
          current_target[toffset + 3 * cont2csize + 25] = current_source[soffset + 253];
          current_target[toffset + 4 * cont2csize + 25] = current_source[soffset + 254];
          current_target[toffset + 5 * cont2csize + 25] = current_source[soffset + 255];
          current_target[toffset + 6 * cont2csize + 25] = current_source[soffset + 256];
          current_target[toffset + 7 * cont2csize + 25] = current_source[soffset + 257];
          current_target[toffset + 8 * cont2csize + 25] = current_source[soffset + 258];
          current_target[toffset + 9 * cont2csize + 25] = current_source[soffset + 259];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 260];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 261];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 262];
          current_target[toffset + 3 * cont2csize + 26] = current_source[soffset + 263];
          current_target[toffset + 4 * cont2csize + 26] = current_source[soffset + 264];
          current_target[toffset + 5 * cont2csize + 26] = current_source[soffset + 265];
          current_target[toffset + 6 * cont2csize + 26] = current_source[soffset + 266];
          current_target[toffset + 7 * cont2csize + 26] = current_source[soffset + 267];
          current_target[toffset + 8 * cont2csize + 26] = current_source[soffset + 268];
          current_target[toffset + 9 * cont2csize + 26] = current_source[soffset + 269];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 270];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 271];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 272];
          current_target[toffset + 3 * cont2csize + 27] = current_source[soffset + 273];
          current_target[toffset + 4 * cont2csize + 27] = current_source[soffset + 274];
          current_target[toffset + 5 * cont2csize + 27] = current_source[soffset + 275];
          current_target[toffset + 6 * cont2csize + 27] = current_source[soffset + 276];
          current_target[toffset + 7 * cont2csize + 27] = current_source[soffset + 277];
          current_target[toffset + 8 * cont2csize + 27] = current_source[soffset + 278];
          current_target[toffset + 9 * cont2csize + 27] = current_source[soffset + 279];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 10 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 10;
          const int soffset = 280 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  10, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 10,  10, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 20,  10, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 30,  10, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 40,  10, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 50,  10, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 60,  10, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+ 70,  10, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+ 80,  10, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+ 90,  10, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+100,  10, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+110,  10, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+120,  10, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+130,  10, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+140,  10, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+150,  10, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+160,  10, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+170,  10, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+180,  10, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+190,  10, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+200,  10, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+210,  10, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+220,  10, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+230,  10, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+240,  10, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+250,  10, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+260,  10, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+270,  10, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_46(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 420;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 420 * (c3 + c3end * c2);
          const int toffset = 15 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 15];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 16];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 17];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 18];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 19];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 20];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 21];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 22];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 23];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 24];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 25];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 26];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 27];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 30];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 31];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 32];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 33];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 34];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 35];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 36];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 37];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 38];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 39];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 40];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 41];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 42];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 43];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 44];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 45];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 46];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 47];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 48];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 49];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 50];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 51];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 52];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 53];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 54];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 55];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 56];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 57];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 58];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 59];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 60];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 61];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 62];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 63];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 64];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 65];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 66];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 67];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 68];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 69];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 70];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 71];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 72];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 73];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 74];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 75];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 76];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 77];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 78];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 79];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 80];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 81];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 82];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 83];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 84];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 85];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 86];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 87];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 88];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 89];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 90];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 91];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 92];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 93];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 94];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 95];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 96];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 97];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 98];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 99];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 100];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 101];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 102];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 103];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 104];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 105];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 106];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 107];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 108];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 109];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 110];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 111];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 112];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 113];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 114];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 115];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 116];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 117];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 118];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 119];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 120];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 121];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 122];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 123];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 124];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 125];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 126];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 127];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 128];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 129];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 130];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 131];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 132];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 133];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 134];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 135];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 136];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 137];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 138];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 139];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 140];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 141];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 142];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 143];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 144];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 145];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 146];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 147];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 148];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 149];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 150];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 151];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 152];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 153];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 154];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 155];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 156];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 157];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 158];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 159];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 160];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 161];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 162];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 163];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 164];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 165];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 166];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 167];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 168];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 169];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 170];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 171];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 172];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 173];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 174];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 175];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 176];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 177];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 178];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 179];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 180];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 181];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 182];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 183];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 184];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 185];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 186];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 187];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 188];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 189];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 190];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 191];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 192];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 193];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 194];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 195];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 196];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 197];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 198];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 199];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 200];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 201];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 202];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 203];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 204];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 205];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 206];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 207];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 208];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 219];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 220];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 221];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 222];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 223];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 224];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 225];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 226];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 227];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 228];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 229];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 230];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 231];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 232];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 233];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 234];
          current_target[toffset + 10 * cont2csize + 15] = current_source[soffset + 235];
          current_target[toffset + 11 * cont2csize + 15] = current_source[soffset + 236];
          current_target[toffset + 12 * cont2csize + 15] = current_source[soffset + 237];
          current_target[toffset + 13 * cont2csize + 15] = current_source[soffset + 238];
          current_target[toffset + 14 * cont2csize + 15] = current_source[soffset + 239];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 240];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 241];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 242];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 243];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 244];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 245];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 246];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 247];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 248];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 249];
          current_target[toffset + 10 * cont2csize + 16] = current_source[soffset + 250];
          current_target[toffset + 11 * cont2csize + 16] = current_source[soffset + 251];
          current_target[toffset + 12 * cont2csize + 16] = current_source[soffset + 252];
          current_target[toffset + 13 * cont2csize + 16] = current_source[soffset + 253];
          current_target[toffset + 14 * cont2csize + 16] = current_source[soffset + 254];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 255];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 256];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 257];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 258];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 259];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 260];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 261];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 262];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 263];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 264];
          current_target[toffset + 10 * cont2csize + 17] = current_source[soffset + 265];
          current_target[toffset + 11 * cont2csize + 17] = current_source[soffset + 266];
          current_target[toffset + 12 * cont2csize + 17] = current_source[soffset + 267];
          current_target[toffset + 13 * cont2csize + 17] = current_source[soffset + 268];
          current_target[toffset + 14 * cont2csize + 17] = current_source[soffset + 269];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 270];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 271];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 272];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 273];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 274];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 275];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 276];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 277];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 278];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 279];
          current_target[toffset + 10 * cont2csize + 18] = current_source[soffset + 280];
          current_target[toffset + 11 * cont2csize + 18] = current_source[soffset + 281];
          current_target[toffset + 12 * cont2csize + 18] = current_source[soffset + 282];
          current_target[toffset + 13 * cont2csize + 18] = current_source[soffset + 283];
          current_target[toffset + 14 * cont2csize + 18] = current_source[soffset + 284];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 285];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 286];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 287];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 288];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 289];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 290];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 291];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 292];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 293];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 294];
          current_target[toffset + 10 * cont2csize + 19] = current_source[soffset + 295];
          current_target[toffset + 11 * cont2csize + 19] = current_source[soffset + 296];
          current_target[toffset + 12 * cont2csize + 19] = current_source[soffset + 297];
          current_target[toffset + 13 * cont2csize + 19] = current_source[soffset + 298];
          current_target[toffset + 14 * cont2csize + 19] = current_source[soffset + 299];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 300];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 301];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 302];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 303];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 304];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 305];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 306];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 307];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 308];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 309];
          current_target[toffset + 10 * cont2csize + 20] = current_source[soffset + 310];
          current_target[toffset + 11 * cont2csize + 20] = current_source[soffset + 311];
          current_target[toffset + 12 * cont2csize + 20] = current_source[soffset + 312];
          current_target[toffset + 13 * cont2csize + 20] = current_source[soffset + 313];
          current_target[toffset + 14 * cont2csize + 20] = current_source[soffset + 314];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 315];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 316];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 317];
          current_target[toffset + 3 * cont2csize + 21] = current_source[soffset + 318];
          current_target[toffset + 4 * cont2csize + 21] = current_source[soffset + 319];
          current_target[toffset + 5 * cont2csize + 21] = current_source[soffset + 320];
          current_target[toffset + 6 * cont2csize + 21] = current_source[soffset + 321];
          current_target[toffset + 7 * cont2csize + 21] = current_source[soffset + 322];
          current_target[toffset + 8 * cont2csize + 21] = current_source[soffset + 323];
          current_target[toffset + 9 * cont2csize + 21] = current_source[soffset + 324];
          current_target[toffset + 10 * cont2csize + 21] = current_source[soffset + 325];
          current_target[toffset + 11 * cont2csize + 21] = current_source[soffset + 326];
          current_target[toffset + 12 * cont2csize + 21] = current_source[soffset + 327];
          current_target[toffset + 13 * cont2csize + 21] = current_source[soffset + 328];
          current_target[toffset + 14 * cont2csize + 21] = current_source[soffset + 329];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 330];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 331];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 332];
          current_target[toffset + 3 * cont2csize + 22] = current_source[soffset + 333];
          current_target[toffset + 4 * cont2csize + 22] = current_source[soffset + 334];
          current_target[toffset + 5 * cont2csize + 22] = current_source[soffset + 335];
          current_target[toffset + 6 * cont2csize + 22] = current_source[soffset + 336];
          current_target[toffset + 7 * cont2csize + 22] = current_source[soffset + 337];
          current_target[toffset + 8 * cont2csize + 22] = current_source[soffset + 338];
          current_target[toffset + 9 * cont2csize + 22] = current_source[soffset + 339];
          current_target[toffset + 10 * cont2csize + 22] = current_source[soffset + 340];
          current_target[toffset + 11 * cont2csize + 22] = current_source[soffset + 341];
          current_target[toffset + 12 * cont2csize + 22] = current_source[soffset + 342];
          current_target[toffset + 13 * cont2csize + 22] = current_source[soffset + 343];
          current_target[toffset + 14 * cont2csize + 22] = current_source[soffset + 344];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 345];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 346];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 347];
          current_target[toffset + 3 * cont2csize + 23] = current_source[soffset + 348];
          current_target[toffset + 4 * cont2csize + 23] = current_source[soffset + 349];
          current_target[toffset + 5 * cont2csize + 23] = current_source[soffset + 350];
          current_target[toffset + 6 * cont2csize + 23] = current_source[soffset + 351];
          current_target[toffset + 7 * cont2csize + 23] = current_source[soffset + 352];
          current_target[toffset + 8 * cont2csize + 23] = current_source[soffset + 353];
          current_target[toffset + 9 * cont2csize + 23] = current_source[soffset + 354];
          current_target[toffset + 10 * cont2csize + 23] = current_source[soffset + 355];
          current_target[toffset + 11 * cont2csize + 23] = current_source[soffset + 356];
          current_target[toffset + 12 * cont2csize + 23] = current_source[soffset + 357];
          current_target[toffset + 13 * cont2csize + 23] = current_source[soffset + 358];
          current_target[toffset + 14 * cont2csize + 23] = current_source[soffset + 359];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 360];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 361];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 362];
          current_target[toffset + 3 * cont2csize + 24] = current_source[soffset + 363];
          current_target[toffset + 4 * cont2csize + 24] = current_source[soffset + 364];
          current_target[toffset + 5 * cont2csize + 24] = current_source[soffset + 365];
          current_target[toffset + 6 * cont2csize + 24] = current_source[soffset + 366];
          current_target[toffset + 7 * cont2csize + 24] = current_source[soffset + 367];
          current_target[toffset + 8 * cont2csize + 24] = current_source[soffset + 368];
          current_target[toffset + 9 * cont2csize + 24] = current_source[soffset + 369];
          current_target[toffset + 10 * cont2csize + 24] = current_source[soffset + 370];
          current_target[toffset + 11 * cont2csize + 24] = current_source[soffset + 371];
          current_target[toffset + 12 * cont2csize + 24] = current_source[soffset + 372];
          current_target[toffset + 13 * cont2csize + 24] = current_source[soffset + 373];
          current_target[toffset + 14 * cont2csize + 24] = current_source[soffset + 374];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 375];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 376];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 377];
          current_target[toffset + 3 * cont2csize + 25] = current_source[soffset + 378];
          current_target[toffset + 4 * cont2csize + 25] = current_source[soffset + 379];
          current_target[toffset + 5 * cont2csize + 25] = current_source[soffset + 380];
          current_target[toffset + 6 * cont2csize + 25] = current_source[soffset + 381];
          current_target[toffset + 7 * cont2csize + 25] = current_source[soffset + 382];
          current_target[toffset + 8 * cont2csize + 25] = current_source[soffset + 383];
          current_target[toffset + 9 * cont2csize + 25] = current_source[soffset + 384];
          current_target[toffset + 10 * cont2csize + 25] = current_source[soffset + 385];
          current_target[toffset + 11 * cont2csize + 25] = current_source[soffset + 386];
          current_target[toffset + 12 * cont2csize + 25] = current_source[soffset + 387];
          current_target[toffset + 13 * cont2csize + 25] = current_source[soffset + 388];
          current_target[toffset + 14 * cont2csize + 25] = current_source[soffset + 389];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 390];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 391];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 392];
          current_target[toffset + 3 * cont2csize + 26] = current_source[soffset + 393];
          current_target[toffset + 4 * cont2csize + 26] = current_source[soffset + 394];
          current_target[toffset + 5 * cont2csize + 26] = current_source[soffset + 395];
          current_target[toffset + 6 * cont2csize + 26] = current_source[soffset + 396];
          current_target[toffset + 7 * cont2csize + 26] = current_source[soffset + 397];
          current_target[toffset + 8 * cont2csize + 26] = current_source[soffset + 398];
          current_target[toffset + 9 * cont2csize + 26] = current_source[soffset + 399];
          current_target[toffset + 10 * cont2csize + 26] = current_source[soffset + 400];
          current_target[toffset + 11 * cont2csize + 26] = current_source[soffset + 401];
          current_target[toffset + 12 * cont2csize + 26] = current_source[soffset + 402];
          current_target[toffset + 13 * cont2csize + 26] = current_source[soffset + 403];
          current_target[toffset + 14 * cont2csize + 26] = current_source[soffset + 404];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 405];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 406];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 407];
          current_target[toffset + 3 * cont2csize + 27] = current_source[soffset + 408];
          current_target[toffset + 4 * cont2csize + 27] = current_source[soffset + 409];
          current_target[toffset + 5 * cont2csize + 27] = current_source[soffset + 410];
          current_target[toffset + 6 * cont2csize + 27] = current_source[soffset + 411];
          current_target[toffset + 7 * cont2csize + 27] = current_source[soffset + 412];
          current_target[toffset + 8 * cont2csize + 27] = current_source[soffset + 413];
          current_target[toffset + 9 * cont2csize + 27] = current_source[soffset + 414];
          current_target[toffset + 10 * cont2csize + 27] = current_source[soffset + 415];
          current_target[toffset + 11 * cont2csize + 27] = current_source[soffset + 416];
          current_target[toffset + 12 * cont2csize + 27] = current_source[soffset + 417];
          current_target[toffset + 13 * cont2csize + 27] = current_source[soffset + 418];
          current_target[toffset + 14 * cont2csize + 27] = current_source[soffset + 419];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 15 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 15;
          const int soffset = 420 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  15, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 15,  15, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 30,  15, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 45,  15, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 60,  15, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+ 75,  15, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+ 90,  15, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+105,  15, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+120,  15, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+135,  15, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+150,  15, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+165,  15, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+180,  15, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+195,  15, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+210,  15, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+225,  15, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+240,  15, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+255,  15, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+270,  15, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+285,  15, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+300,  15, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+315,  15, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+330,  15, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+345,  15, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+360,  15, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+375,  15, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+390,  15, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+405,  15, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_56(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 588;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 588 * (c3 + c3end * c2);
          const int toffset = 21 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 15 * cont2csize + 0] = current_source[soffset + 15];
          current_target[toffset + 16 * cont2csize + 0] = current_source[soffset + 16];
          current_target[toffset + 17 * cont2csize + 0] = current_source[soffset + 17];
          current_target[toffset + 18 * cont2csize + 0] = current_source[soffset + 18];
          current_target[toffset + 19 * cont2csize + 0] = current_source[soffset + 19];
          current_target[toffset + 20 * cont2csize + 0] = current_source[soffset + 20];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 21];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 22];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 23];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 24];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 25];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 26];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 27];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 30];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 31];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 32];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 33];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 34];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 35];
          current_target[toffset + 15 * cont2csize + 1] = current_source[soffset + 36];
          current_target[toffset + 16 * cont2csize + 1] = current_source[soffset + 37];
          current_target[toffset + 17 * cont2csize + 1] = current_source[soffset + 38];
          current_target[toffset + 18 * cont2csize + 1] = current_source[soffset + 39];
          current_target[toffset + 19 * cont2csize + 1] = current_source[soffset + 40];
          current_target[toffset + 20 * cont2csize + 1] = current_source[soffset + 41];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 42];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 43];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 44];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 45];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 46];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 47];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 48];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 49];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 50];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 51];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 52];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 53];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 54];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 55];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 56];
          current_target[toffset + 15 * cont2csize + 2] = current_source[soffset + 57];
          current_target[toffset + 16 * cont2csize + 2] = current_source[soffset + 58];
          current_target[toffset + 17 * cont2csize + 2] = current_source[soffset + 59];
          current_target[toffset + 18 * cont2csize + 2] = current_source[soffset + 60];
          current_target[toffset + 19 * cont2csize + 2] = current_source[soffset + 61];
          current_target[toffset + 20 * cont2csize + 2] = current_source[soffset + 62];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 63];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 64];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 65];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 66];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 67];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 68];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 69];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 70];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 71];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 72];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 73];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 74];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 75];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 76];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 77];
          current_target[toffset + 15 * cont2csize + 3] = current_source[soffset + 78];
          current_target[toffset + 16 * cont2csize + 3] = current_source[soffset + 79];
          current_target[toffset + 17 * cont2csize + 3] = current_source[soffset + 80];
          current_target[toffset + 18 * cont2csize + 3] = current_source[soffset + 81];
          current_target[toffset + 19 * cont2csize + 3] = current_source[soffset + 82];
          current_target[toffset + 20 * cont2csize + 3] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 89];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 90];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 91];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 92];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 93];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 94];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 95];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 96];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 97];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 98];
          current_target[toffset + 15 * cont2csize + 4] = current_source[soffset + 99];
          current_target[toffset + 16 * cont2csize + 4] = current_source[soffset + 100];
          current_target[toffset + 17 * cont2csize + 4] = current_source[soffset + 101];
          current_target[toffset + 18 * cont2csize + 4] = current_source[soffset + 102];
          current_target[toffset + 19 * cont2csize + 4] = current_source[soffset + 103];
          current_target[toffset + 20 * cont2csize + 4] = current_source[soffset + 104];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 105];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 106];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 107];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 108];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 109];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 110];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 111];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 112];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 113];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 114];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 115];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 116];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 117];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 118];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 119];
          current_target[toffset + 15 * cont2csize + 5] = current_source[soffset + 120];
          current_target[toffset + 16 * cont2csize + 5] = current_source[soffset + 121];
          current_target[toffset + 17 * cont2csize + 5] = current_source[soffset + 122];
          current_target[toffset + 18 * cont2csize + 5] = current_source[soffset + 123];
          current_target[toffset + 19 * cont2csize + 5] = current_source[soffset + 124];
          current_target[toffset + 20 * cont2csize + 5] = current_source[soffset + 125];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 126];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 127];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 128];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 129];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 130];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 131];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 132];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 133];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 134];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 135];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 136];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 137];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 138];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 139];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 140];
          current_target[toffset + 15 * cont2csize + 6] = current_source[soffset + 141];
          current_target[toffset + 16 * cont2csize + 6] = current_source[soffset + 142];
          current_target[toffset + 17 * cont2csize + 6] = current_source[soffset + 143];
          current_target[toffset + 18 * cont2csize + 6] = current_source[soffset + 144];
          current_target[toffset + 19 * cont2csize + 6] = current_source[soffset + 145];
          current_target[toffset + 20 * cont2csize + 6] = current_source[soffset + 146];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 147];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 148];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 149];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 150];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 151];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 152];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 153];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 154];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 155];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 156];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 157];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 158];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 159];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 160];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 161];
          current_target[toffset + 15 * cont2csize + 7] = current_source[soffset + 162];
          current_target[toffset + 16 * cont2csize + 7] = current_source[soffset + 163];
          current_target[toffset + 17 * cont2csize + 7] = current_source[soffset + 164];
          current_target[toffset + 18 * cont2csize + 7] = current_source[soffset + 165];
          current_target[toffset + 19 * cont2csize + 7] = current_source[soffset + 166];
          current_target[toffset + 20 * cont2csize + 7] = current_source[soffset + 167];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 168];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 169];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 170];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 171];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 172];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 173];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 174];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 175];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 176];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 177];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 178];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 179];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 180];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 181];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 182];
          current_target[toffset + 15 * cont2csize + 8] = current_source[soffset + 183];
          current_target[toffset + 16 * cont2csize + 8] = current_source[soffset + 184];
          current_target[toffset + 17 * cont2csize + 8] = current_source[soffset + 185];
          current_target[toffset + 18 * cont2csize + 8] = current_source[soffset + 186];
          current_target[toffset + 19 * cont2csize + 8] = current_source[soffset + 187];
          current_target[toffset + 20 * cont2csize + 8] = current_source[soffset + 188];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 189];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 190];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 191];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 192];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 193];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 194];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 195];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 196];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 197];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 198];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 199];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 200];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 201];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 202];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 203];
          current_target[toffset + 15 * cont2csize + 9] = current_source[soffset + 204];
          current_target[toffset + 16 * cont2csize + 9] = current_source[soffset + 205];
          current_target[toffset + 17 * cont2csize + 9] = current_source[soffset + 206];
          current_target[toffset + 18 * cont2csize + 9] = current_source[soffset + 207];
          current_target[toffset + 19 * cont2csize + 9] = current_source[soffset + 208];
          current_target[toffset + 20 * cont2csize + 9] = current_source[soffset + 209];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 210];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 211];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 212];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 213];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 214];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 215];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 216];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 217];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 218];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 219];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 220];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 221];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 222];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 223];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 224];
          current_target[toffset + 15 * cont2csize + 10] = current_source[soffset + 225];
          current_target[toffset + 16 * cont2csize + 10] = current_source[soffset + 226];
          current_target[toffset + 17 * cont2csize + 10] = current_source[soffset + 227];
          current_target[toffset + 18 * cont2csize + 10] = current_source[soffset + 228];
          current_target[toffset + 19 * cont2csize + 10] = current_source[soffset + 229];
          current_target[toffset + 20 * cont2csize + 10] = current_source[soffset + 230];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 231];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 232];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 233];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 234];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 235];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 236];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 237];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 238];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 239];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 240];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 241];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 242];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 243];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 244];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 245];
          current_target[toffset + 15 * cont2csize + 11] = current_source[soffset + 246];
          current_target[toffset + 16 * cont2csize + 11] = current_source[soffset + 247];
          current_target[toffset + 17 * cont2csize + 11] = current_source[soffset + 248];
          current_target[toffset + 18 * cont2csize + 11] = current_source[soffset + 249];
          current_target[toffset + 19 * cont2csize + 11] = current_source[soffset + 250];
          current_target[toffset + 20 * cont2csize + 11] = current_source[soffset + 251];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 252];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 253];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 254];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 255];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 256];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 257];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 258];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 259];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 260];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 261];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 262];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 263];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 264];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 265];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 266];
          current_target[toffset + 15 * cont2csize + 12] = current_source[soffset + 267];
          current_target[toffset + 16 * cont2csize + 12] = current_source[soffset + 268];
          current_target[toffset + 17 * cont2csize + 12] = current_source[soffset + 269];
          current_target[toffset + 18 * cont2csize + 12] = current_source[soffset + 270];
          current_target[toffset + 19 * cont2csize + 12] = current_source[soffset + 271];
          current_target[toffset + 20 * cont2csize + 12] = current_source[soffset + 272];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 273];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 274];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 275];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 276];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 277];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 278];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 279];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 280];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 281];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 282];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 283];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 284];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 285];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 286];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 287];
          current_target[toffset + 15 * cont2csize + 13] = current_source[soffset + 288];
          current_target[toffset + 16 * cont2csize + 13] = current_source[soffset + 289];
          current_target[toffset + 17 * cont2csize + 13] = current_source[soffset + 290];
          current_target[toffset + 18 * cont2csize + 13] = current_source[soffset + 291];
          current_target[toffset + 19 * cont2csize + 13] = current_source[soffset + 292];
          current_target[toffset + 20 * cont2csize + 13] = current_source[soffset + 293];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 294];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 295];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 296];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 297];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 298];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 299];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 300];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 301];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 302];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 303];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 304];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 305];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 306];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 307];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 308];
          current_target[toffset + 15 * cont2csize + 14] = current_source[soffset + 309];
          current_target[toffset + 16 * cont2csize + 14] = current_source[soffset + 310];
          current_target[toffset + 17 * cont2csize + 14] = current_source[soffset + 311];
          current_target[toffset + 18 * cont2csize + 14] = current_source[soffset + 312];
          current_target[toffset + 19 * cont2csize + 14] = current_source[soffset + 313];
          current_target[toffset + 20 * cont2csize + 14] = current_source[soffset + 314];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 315];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 316];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 317];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 318];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 319];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 320];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 321];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 322];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 323];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 324];
          current_target[toffset + 10 * cont2csize + 15] = current_source[soffset + 325];
          current_target[toffset + 11 * cont2csize + 15] = current_source[soffset + 326];
          current_target[toffset + 12 * cont2csize + 15] = current_source[soffset + 327];
          current_target[toffset + 13 * cont2csize + 15] = current_source[soffset + 328];
          current_target[toffset + 14 * cont2csize + 15] = current_source[soffset + 329];
          current_target[toffset + 15 * cont2csize + 15] = current_source[soffset + 330];
          current_target[toffset + 16 * cont2csize + 15] = current_source[soffset + 331];
          current_target[toffset + 17 * cont2csize + 15] = current_source[soffset + 332];
          current_target[toffset + 18 * cont2csize + 15] = current_source[soffset + 333];
          current_target[toffset + 19 * cont2csize + 15] = current_source[soffset + 334];
          current_target[toffset + 20 * cont2csize + 15] = current_source[soffset + 335];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 336];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 337];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 338];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 339];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 340];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 341];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 342];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 343];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 344];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 345];
          current_target[toffset + 10 * cont2csize + 16] = current_source[soffset + 346];
          current_target[toffset + 11 * cont2csize + 16] = current_source[soffset + 347];
          current_target[toffset + 12 * cont2csize + 16] = current_source[soffset + 348];
          current_target[toffset + 13 * cont2csize + 16] = current_source[soffset + 349];
          current_target[toffset + 14 * cont2csize + 16] = current_source[soffset + 350];
          current_target[toffset + 15 * cont2csize + 16] = current_source[soffset + 351];
          current_target[toffset + 16 * cont2csize + 16] = current_source[soffset + 352];
          current_target[toffset + 17 * cont2csize + 16] = current_source[soffset + 353];
          current_target[toffset + 18 * cont2csize + 16] = current_source[soffset + 354];
          current_target[toffset + 19 * cont2csize + 16] = current_source[soffset + 355];
          current_target[toffset + 20 * cont2csize + 16] = current_source[soffset + 356];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 357];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 358];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 359];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 360];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 361];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 362];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 363];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 364];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 365];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 366];
          current_target[toffset + 10 * cont2csize + 17] = current_source[soffset + 367];
          current_target[toffset + 11 * cont2csize + 17] = current_source[soffset + 368];
          current_target[toffset + 12 * cont2csize + 17] = current_source[soffset + 369];
          current_target[toffset + 13 * cont2csize + 17] = current_source[soffset + 370];
          current_target[toffset + 14 * cont2csize + 17] = current_source[soffset + 371];
          current_target[toffset + 15 * cont2csize + 17] = current_source[soffset + 372];
          current_target[toffset + 16 * cont2csize + 17] = current_source[soffset + 373];
          current_target[toffset + 17 * cont2csize + 17] = current_source[soffset + 374];
          current_target[toffset + 18 * cont2csize + 17] = current_source[soffset + 375];
          current_target[toffset + 19 * cont2csize + 17] = current_source[soffset + 376];
          current_target[toffset + 20 * cont2csize + 17] = current_source[soffset + 377];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 378];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 379];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 380];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 381];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 382];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 383];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 384];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 385];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 386];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 387];
          current_target[toffset + 10 * cont2csize + 18] = current_source[soffset + 388];
          current_target[toffset + 11 * cont2csize + 18] = current_source[soffset + 389];
          current_target[toffset + 12 * cont2csize + 18] = current_source[soffset + 390];
          current_target[toffset + 13 * cont2csize + 18] = current_source[soffset + 391];
          current_target[toffset + 14 * cont2csize + 18] = current_source[soffset + 392];
          current_target[toffset + 15 * cont2csize + 18] = current_source[soffset + 393];
          current_target[toffset + 16 * cont2csize + 18] = current_source[soffset + 394];
          current_target[toffset + 17 * cont2csize + 18] = current_source[soffset + 395];
          current_target[toffset + 18 * cont2csize + 18] = current_source[soffset + 396];
          current_target[toffset + 19 * cont2csize + 18] = current_source[soffset + 397];
          current_target[toffset + 20 * cont2csize + 18] = current_source[soffset + 398];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 399];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 400];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 401];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 402];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 403];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 404];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 405];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 406];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 407];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 408];
          current_target[toffset + 10 * cont2csize + 19] = current_source[soffset + 409];
          current_target[toffset + 11 * cont2csize + 19] = current_source[soffset + 410];
          current_target[toffset + 12 * cont2csize + 19] = current_source[soffset + 411];
          current_target[toffset + 13 * cont2csize + 19] = current_source[soffset + 412];
          current_target[toffset + 14 * cont2csize + 19] = current_source[soffset + 413];
          current_target[toffset + 15 * cont2csize + 19] = current_source[soffset + 414];
          current_target[toffset + 16 * cont2csize + 19] = current_source[soffset + 415];
          current_target[toffset + 17 * cont2csize + 19] = current_source[soffset + 416];
          current_target[toffset + 18 * cont2csize + 19] = current_source[soffset + 417];
          current_target[toffset + 19 * cont2csize + 19] = current_source[soffset + 418];
          current_target[toffset + 20 * cont2csize + 19] = current_source[soffset + 419];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 420];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 421];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 422];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 423];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 424];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 425];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 426];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 427];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 428];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 429];
          current_target[toffset + 10 * cont2csize + 20] = current_source[soffset + 430];
          current_target[toffset + 11 * cont2csize + 20] = current_source[soffset + 431];
          current_target[toffset + 12 * cont2csize + 20] = current_source[soffset + 432];
          current_target[toffset + 13 * cont2csize + 20] = current_source[soffset + 433];
          current_target[toffset + 14 * cont2csize + 20] = current_source[soffset + 434];
          current_target[toffset + 15 * cont2csize + 20] = current_source[soffset + 435];
          current_target[toffset + 16 * cont2csize + 20] = current_source[soffset + 436];
          current_target[toffset + 17 * cont2csize + 20] = current_source[soffset + 437];
          current_target[toffset + 18 * cont2csize + 20] = current_source[soffset + 438];
          current_target[toffset + 19 * cont2csize + 20] = current_source[soffset + 439];
          current_target[toffset + 20 * cont2csize + 20] = current_source[soffset + 440];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 441];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 442];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 443];
          current_target[toffset + 3 * cont2csize + 21] = current_source[soffset + 444];
          current_target[toffset + 4 * cont2csize + 21] = current_source[soffset + 445];
          current_target[toffset + 5 * cont2csize + 21] = current_source[soffset + 446];
          current_target[toffset + 6 * cont2csize + 21] = current_source[soffset + 447];
          current_target[toffset + 7 * cont2csize + 21] = current_source[soffset + 448];
          current_target[toffset + 8 * cont2csize + 21] = current_source[soffset + 449];
          current_target[toffset + 9 * cont2csize + 21] = current_source[soffset + 450];
          current_target[toffset + 10 * cont2csize + 21] = current_source[soffset + 451];
          current_target[toffset + 11 * cont2csize + 21] = current_source[soffset + 452];
          current_target[toffset + 12 * cont2csize + 21] = current_source[soffset + 453];
          current_target[toffset + 13 * cont2csize + 21] = current_source[soffset + 454];
          current_target[toffset + 14 * cont2csize + 21] = current_source[soffset + 455];
          current_target[toffset + 15 * cont2csize + 21] = current_source[soffset + 456];
          current_target[toffset + 16 * cont2csize + 21] = current_source[soffset + 457];
          current_target[toffset + 17 * cont2csize + 21] = current_source[soffset + 458];
          current_target[toffset + 18 * cont2csize + 21] = current_source[soffset + 459];
          current_target[toffset + 19 * cont2csize + 21] = current_source[soffset + 460];
          current_target[toffset + 20 * cont2csize + 21] = current_source[soffset + 461];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 462];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 463];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 464];
          current_target[toffset + 3 * cont2csize + 22] = current_source[soffset + 465];
          current_target[toffset + 4 * cont2csize + 22] = current_source[soffset + 466];
          current_target[toffset + 5 * cont2csize + 22] = current_source[soffset + 467];
          current_target[toffset + 6 * cont2csize + 22] = current_source[soffset + 468];
          current_target[toffset + 7 * cont2csize + 22] = current_source[soffset + 469];
          current_target[toffset + 8 * cont2csize + 22] = current_source[soffset + 470];
          current_target[toffset + 9 * cont2csize + 22] = current_source[soffset + 471];
          current_target[toffset + 10 * cont2csize + 22] = current_source[soffset + 472];
          current_target[toffset + 11 * cont2csize + 22] = current_source[soffset + 473];
          current_target[toffset + 12 * cont2csize + 22] = current_source[soffset + 474];
          current_target[toffset + 13 * cont2csize + 22] = current_source[soffset + 475];
          current_target[toffset + 14 * cont2csize + 22] = current_source[soffset + 476];
          current_target[toffset + 15 * cont2csize + 22] = current_source[soffset + 477];
          current_target[toffset + 16 * cont2csize + 22] = current_source[soffset + 478];
          current_target[toffset + 17 * cont2csize + 22] = current_source[soffset + 479];
          current_target[toffset + 18 * cont2csize + 22] = current_source[soffset + 480];
          current_target[toffset + 19 * cont2csize + 22] = current_source[soffset + 481];
          current_target[toffset + 20 * cont2csize + 22] = current_source[soffset + 482];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 483];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 484];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 485];
          current_target[toffset + 3 * cont2csize + 23] = current_source[soffset + 486];
          current_target[toffset + 4 * cont2csize + 23] = current_source[soffset + 487];
          current_target[toffset + 5 * cont2csize + 23] = current_source[soffset + 488];
          current_target[toffset + 6 * cont2csize + 23] = current_source[soffset + 489];
          current_target[toffset + 7 * cont2csize + 23] = current_source[soffset + 490];
          current_target[toffset + 8 * cont2csize + 23] = current_source[soffset + 491];
          current_target[toffset + 9 * cont2csize + 23] = current_source[soffset + 492];
          current_target[toffset + 10 * cont2csize + 23] = current_source[soffset + 493];
          current_target[toffset + 11 * cont2csize + 23] = current_source[soffset + 494];
          current_target[toffset + 12 * cont2csize + 23] = current_source[soffset + 495];
          current_target[toffset + 13 * cont2csize + 23] = current_source[soffset + 496];
          current_target[toffset + 14 * cont2csize + 23] = current_source[soffset + 497];
          current_target[toffset + 15 * cont2csize + 23] = current_source[soffset + 498];
          current_target[toffset + 16 * cont2csize + 23] = current_source[soffset + 499];
          current_target[toffset + 17 * cont2csize + 23] = current_source[soffset + 500];
          current_target[toffset + 18 * cont2csize + 23] = current_source[soffset + 501];
          current_target[toffset + 19 * cont2csize + 23] = current_source[soffset + 502];
          current_target[toffset + 20 * cont2csize + 23] = current_source[soffset + 503];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 504];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 505];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 506];
          current_target[toffset + 3 * cont2csize + 24] = current_source[soffset + 507];
          current_target[toffset + 4 * cont2csize + 24] = current_source[soffset + 508];
          current_target[toffset + 5 * cont2csize + 24] = current_source[soffset + 509];
          current_target[toffset + 6 * cont2csize + 24] = current_source[soffset + 510];
          current_target[toffset + 7 * cont2csize + 24] = current_source[soffset + 511];
          current_target[toffset + 8 * cont2csize + 24] = current_source[soffset + 512];
          current_target[toffset + 9 * cont2csize + 24] = current_source[soffset + 513];
          current_target[toffset + 10 * cont2csize + 24] = current_source[soffset + 514];
          current_target[toffset + 11 * cont2csize + 24] = current_source[soffset + 515];
          current_target[toffset + 12 * cont2csize + 24] = current_source[soffset + 516];
          current_target[toffset + 13 * cont2csize + 24] = current_source[soffset + 517];
          current_target[toffset + 14 * cont2csize + 24] = current_source[soffset + 518];
          current_target[toffset + 15 * cont2csize + 24] = current_source[soffset + 519];
          current_target[toffset + 16 * cont2csize + 24] = current_source[soffset + 520];
          current_target[toffset + 17 * cont2csize + 24] = current_source[soffset + 521];
          current_target[toffset + 18 * cont2csize + 24] = current_source[soffset + 522];
          current_target[toffset + 19 * cont2csize + 24] = current_source[soffset + 523];
          current_target[toffset + 20 * cont2csize + 24] = current_source[soffset + 524];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 525];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 526];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 527];
          current_target[toffset + 3 * cont2csize + 25] = current_source[soffset + 528];
          current_target[toffset + 4 * cont2csize + 25] = current_source[soffset + 529];
          current_target[toffset + 5 * cont2csize + 25] = current_source[soffset + 530];
          current_target[toffset + 6 * cont2csize + 25] = current_source[soffset + 531];
          current_target[toffset + 7 * cont2csize + 25] = current_source[soffset + 532];
          current_target[toffset + 8 * cont2csize + 25] = current_source[soffset + 533];
          current_target[toffset + 9 * cont2csize + 25] = current_source[soffset + 534];
          current_target[toffset + 10 * cont2csize + 25] = current_source[soffset + 535];
          current_target[toffset + 11 * cont2csize + 25] = current_source[soffset + 536];
          current_target[toffset + 12 * cont2csize + 25] = current_source[soffset + 537];
          current_target[toffset + 13 * cont2csize + 25] = current_source[soffset + 538];
          current_target[toffset + 14 * cont2csize + 25] = current_source[soffset + 539];
          current_target[toffset + 15 * cont2csize + 25] = current_source[soffset + 540];
          current_target[toffset + 16 * cont2csize + 25] = current_source[soffset + 541];
          current_target[toffset + 17 * cont2csize + 25] = current_source[soffset + 542];
          current_target[toffset + 18 * cont2csize + 25] = current_source[soffset + 543];
          current_target[toffset + 19 * cont2csize + 25] = current_source[soffset + 544];
          current_target[toffset + 20 * cont2csize + 25] = current_source[soffset + 545];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 546];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 547];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 548];
          current_target[toffset + 3 * cont2csize + 26] = current_source[soffset + 549];
          current_target[toffset + 4 * cont2csize + 26] = current_source[soffset + 550];
          current_target[toffset + 5 * cont2csize + 26] = current_source[soffset + 551];
          current_target[toffset + 6 * cont2csize + 26] = current_source[soffset + 552];
          current_target[toffset + 7 * cont2csize + 26] = current_source[soffset + 553];
          current_target[toffset + 8 * cont2csize + 26] = current_source[soffset + 554];
          current_target[toffset + 9 * cont2csize + 26] = current_source[soffset + 555];
          current_target[toffset + 10 * cont2csize + 26] = current_source[soffset + 556];
          current_target[toffset + 11 * cont2csize + 26] = current_source[soffset + 557];
          current_target[toffset + 12 * cont2csize + 26] = current_source[soffset + 558];
          current_target[toffset + 13 * cont2csize + 26] = current_source[soffset + 559];
          current_target[toffset + 14 * cont2csize + 26] = current_source[soffset + 560];
          current_target[toffset + 15 * cont2csize + 26] = current_source[soffset + 561];
          current_target[toffset + 16 * cont2csize + 26] = current_source[soffset + 562];
          current_target[toffset + 17 * cont2csize + 26] = current_source[soffset + 563];
          current_target[toffset + 18 * cont2csize + 26] = current_source[soffset + 564];
          current_target[toffset + 19 * cont2csize + 26] = current_source[soffset + 565];
          current_target[toffset + 20 * cont2csize + 26] = current_source[soffset + 566];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 567];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 568];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 569];
          current_target[toffset + 3 * cont2csize + 27] = current_source[soffset + 570];
          current_target[toffset + 4 * cont2csize + 27] = current_source[soffset + 571];
          current_target[toffset + 5 * cont2csize + 27] = current_source[soffset + 572];
          current_target[toffset + 6 * cont2csize + 27] = current_source[soffset + 573];
          current_target[toffset + 7 * cont2csize + 27] = current_source[soffset + 574];
          current_target[toffset + 8 * cont2csize + 27] = current_source[soffset + 575];
          current_target[toffset + 9 * cont2csize + 27] = current_source[soffset + 576];
          current_target[toffset + 10 * cont2csize + 27] = current_source[soffset + 577];
          current_target[toffset + 11 * cont2csize + 27] = current_source[soffset + 578];
          current_target[toffset + 12 * cont2csize + 27] = current_source[soffset + 579];
          current_target[toffset + 13 * cont2csize + 27] = current_source[soffset + 580];
          current_target[toffset + 14 * cont2csize + 27] = current_source[soffset + 581];
          current_target[toffset + 15 * cont2csize + 27] = current_source[soffset + 582];
          current_target[toffset + 16 * cont2csize + 27] = current_source[soffset + 583];
          current_target[toffset + 17 * cont2csize + 27] = current_source[soffset + 584];
          current_target[toffset + 18 * cont2csize + 27] = current_source[soffset + 585];
          current_target[toffset + 19 * cont2csize + 27] = current_source[soffset + 586];
          current_target[toffset + 20 * cont2csize + 27] = current_source[soffset + 587];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 21 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 21;
          const int soffset = 588 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  21, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 21,  21, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 42,  21, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 63,  21, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+ 84,  21, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+105,  21, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+126,  21, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+147,  21, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+168,  21, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+189,  21, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+210,  21, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+231,  21, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+252,  21, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+273,  21, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+294,  21, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+315,  21, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+336,  21, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+357,  21, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+378,  21, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+399,  21, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+420,  21, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+441,  21, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+462,  21, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+483,  21, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+504,  21, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+525,  21, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+546,  21, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+567,  21, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}


void CSortList::sort_indices_66(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {
  const int innerloopsize = c2end * c3end * 784;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = 28 * c2end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        const int c2x2end = c2 * 28;
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int soffset = 784 * (c3 + c3end * c2);
          const int toffset = 28 * c3 * cont2csize + c2x2end;
          current_target[toffset + 0 * cont2csize + 0] = current_source[soffset + 0];
          current_target[toffset + 1 * cont2csize + 0] = current_source[soffset + 1];
          current_target[toffset + 2 * cont2csize + 0] = current_source[soffset + 2];
          current_target[toffset + 3 * cont2csize + 0] = current_source[soffset + 3];
          current_target[toffset + 4 * cont2csize + 0] = current_source[soffset + 4];
          current_target[toffset + 5 * cont2csize + 0] = current_source[soffset + 5];
          current_target[toffset + 6 * cont2csize + 0] = current_source[soffset + 6];
          current_target[toffset + 7 * cont2csize + 0] = current_source[soffset + 7];
          current_target[toffset + 8 * cont2csize + 0] = current_source[soffset + 8];
          current_target[toffset + 9 * cont2csize + 0] = current_source[soffset + 9];
          current_target[toffset + 10 * cont2csize + 0] = current_source[soffset + 10];
          current_target[toffset + 11 * cont2csize + 0] = current_source[soffset + 11];
          current_target[toffset + 12 * cont2csize + 0] = current_source[soffset + 12];
          current_target[toffset + 13 * cont2csize + 0] = current_source[soffset + 13];
          current_target[toffset + 14 * cont2csize + 0] = current_source[soffset + 14];
          current_target[toffset + 15 * cont2csize + 0] = current_source[soffset + 15];
          current_target[toffset + 16 * cont2csize + 0] = current_source[soffset + 16];
          current_target[toffset + 17 * cont2csize + 0] = current_source[soffset + 17];
          current_target[toffset + 18 * cont2csize + 0] = current_source[soffset + 18];
          current_target[toffset + 19 * cont2csize + 0] = current_source[soffset + 19];
          current_target[toffset + 20 * cont2csize + 0] = current_source[soffset + 20];
          current_target[toffset + 21 * cont2csize + 0] = current_source[soffset + 21];
          current_target[toffset + 22 * cont2csize + 0] = current_source[soffset + 22];
          current_target[toffset + 23 * cont2csize + 0] = current_source[soffset + 23];
          current_target[toffset + 24 * cont2csize + 0] = current_source[soffset + 24];
          current_target[toffset + 25 * cont2csize + 0] = current_source[soffset + 25];
          current_target[toffset + 26 * cont2csize + 0] = current_source[soffset + 26];
          current_target[toffset + 27 * cont2csize + 0] = current_source[soffset + 27];
          current_target[toffset + 0 * cont2csize + 1] = current_source[soffset + 28];
          current_target[toffset + 1 * cont2csize + 1] = current_source[soffset + 29];
          current_target[toffset + 2 * cont2csize + 1] = current_source[soffset + 30];
          current_target[toffset + 3 * cont2csize + 1] = current_source[soffset + 31];
          current_target[toffset + 4 * cont2csize + 1] = current_source[soffset + 32];
          current_target[toffset + 5 * cont2csize + 1] = current_source[soffset + 33];
          current_target[toffset + 6 * cont2csize + 1] = current_source[soffset + 34];
          current_target[toffset + 7 * cont2csize + 1] = current_source[soffset + 35];
          current_target[toffset + 8 * cont2csize + 1] = current_source[soffset + 36];
          current_target[toffset + 9 * cont2csize + 1] = current_source[soffset + 37];
          current_target[toffset + 10 * cont2csize + 1] = current_source[soffset + 38];
          current_target[toffset + 11 * cont2csize + 1] = current_source[soffset + 39];
          current_target[toffset + 12 * cont2csize + 1] = current_source[soffset + 40];
          current_target[toffset + 13 * cont2csize + 1] = current_source[soffset + 41];
          current_target[toffset + 14 * cont2csize + 1] = current_source[soffset + 42];
          current_target[toffset + 15 * cont2csize + 1] = current_source[soffset + 43];
          current_target[toffset + 16 * cont2csize + 1] = current_source[soffset + 44];
          current_target[toffset + 17 * cont2csize + 1] = current_source[soffset + 45];
          current_target[toffset + 18 * cont2csize + 1] = current_source[soffset + 46];
          current_target[toffset + 19 * cont2csize + 1] = current_source[soffset + 47];
          current_target[toffset + 20 * cont2csize + 1] = current_source[soffset + 48];
          current_target[toffset + 21 * cont2csize + 1] = current_source[soffset + 49];
          current_target[toffset + 22 * cont2csize + 1] = current_source[soffset + 50];
          current_target[toffset + 23 * cont2csize + 1] = current_source[soffset + 51];
          current_target[toffset + 24 * cont2csize + 1] = current_source[soffset + 52];
          current_target[toffset + 25 * cont2csize + 1] = current_source[soffset + 53];
          current_target[toffset + 26 * cont2csize + 1] = current_source[soffset + 54];
          current_target[toffset + 27 * cont2csize + 1] = current_source[soffset + 55];
          current_target[toffset + 0 * cont2csize + 2] = current_source[soffset + 56];
          current_target[toffset + 1 * cont2csize + 2] = current_source[soffset + 57];
          current_target[toffset + 2 * cont2csize + 2] = current_source[soffset + 58];
          current_target[toffset + 3 * cont2csize + 2] = current_source[soffset + 59];
          current_target[toffset + 4 * cont2csize + 2] = current_source[soffset + 60];
          current_target[toffset + 5 * cont2csize + 2] = current_source[soffset + 61];
          current_target[toffset + 6 * cont2csize + 2] = current_source[soffset + 62];
          current_target[toffset + 7 * cont2csize + 2] = current_source[soffset + 63];
          current_target[toffset + 8 * cont2csize + 2] = current_source[soffset + 64];
          current_target[toffset + 9 * cont2csize + 2] = current_source[soffset + 65];
          current_target[toffset + 10 * cont2csize + 2] = current_source[soffset + 66];
          current_target[toffset + 11 * cont2csize + 2] = current_source[soffset + 67];
          current_target[toffset + 12 * cont2csize + 2] = current_source[soffset + 68];
          current_target[toffset + 13 * cont2csize + 2] = current_source[soffset + 69];
          current_target[toffset + 14 * cont2csize + 2] = current_source[soffset + 70];
          current_target[toffset + 15 * cont2csize + 2] = current_source[soffset + 71];
          current_target[toffset + 16 * cont2csize + 2] = current_source[soffset + 72];
          current_target[toffset + 17 * cont2csize + 2] = current_source[soffset + 73];
          current_target[toffset + 18 * cont2csize + 2] = current_source[soffset + 74];
          current_target[toffset + 19 * cont2csize + 2] = current_source[soffset + 75];
          current_target[toffset + 20 * cont2csize + 2] = current_source[soffset + 76];
          current_target[toffset + 21 * cont2csize + 2] = current_source[soffset + 77];
          current_target[toffset + 22 * cont2csize + 2] = current_source[soffset + 78];
          current_target[toffset + 23 * cont2csize + 2] = current_source[soffset + 79];
          current_target[toffset + 24 * cont2csize + 2] = current_source[soffset + 80];
          current_target[toffset + 25 * cont2csize + 2] = current_source[soffset + 81];
          current_target[toffset + 26 * cont2csize + 2] = current_source[soffset + 82];
          current_target[toffset + 27 * cont2csize + 2] = current_source[soffset + 83];
          current_target[toffset + 0 * cont2csize + 3] = current_source[soffset + 84];
          current_target[toffset + 1 * cont2csize + 3] = current_source[soffset + 85];
          current_target[toffset + 2 * cont2csize + 3] = current_source[soffset + 86];
          current_target[toffset + 3 * cont2csize + 3] = current_source[soffset + 87];
          current_target[toffset + 4 * cont2csize + 3] = current_source[soffset + 88];
          current_target[toffset + 5 * cont2csize + 3] = current_source[soffset + 89];
          current_target[toffset + 6 * cont2csize + 3] = current_source[soffset + 90];
          current_target[toffset + 7 * cont2csize + 3] = current_source[soffset + 91];
          current_target[toffset + 8 * cont2csize + 3] = current_source[soffset + 92];
          current_target[toffset + 9 * cont2csize + 3] = current_source[soffset + 93];
          current_target[toffset + 10 * cont2csize + 3] = current_source[soffset + 94];
          current_target[toffset + 11 * cont2csize + 3] = current_source[soffset + 95];
          current_target[toffset + 12 * cont2csize + 3] = current_source[soffset + 96];
          current_target[toffset + 13 * cont2csize + 3] = current_source[soffset + 97];
          current_target[toffset + 14 * cont2csize + 3] = current_source[soffset + 98];
          current_target[toffset + 15 * cont2csize + 3] = current_source[soffset + 99];
          current_target[toffset + 16 * cont2csize + 3] = current_source[soffset + 100];
          current_target[toffset + 17 * cont2csize + 3] = current_source[soffset + 101];
          current_target[toffset + 18 * cont2csize + 3] = current_source[soffset + 102];
          current_target[toffset + 19 * cont2csize + 3] = current_source[soffset + 103];
          current_target[toffset + 20 * cont2csize + 3] = current_source[soffset + 104];
          current_target[toffset + 21 * cont2csize + 3] = current_source[soffset + 105];
          current_target[toffset + 22 * cont2csize + 3] = current_source[soffset + 106];
          current_target[toffset + 23 * cont2csize + 3] = current_source[soffset + 107];
          current_target[toffset + 24 * cont2csize + 3] = current_source[soffset + 108];
          current_target[toffset + 25 * cont2csize + 3] = current_source[soffset + 109];
          current_target[toffset + 26 * cont2csize + 3] = current_source[soffset + 110];
          current_target[toffset + 27 * cont2csize + 3] = current_source[soffset + 111];
          current_target[toffset + 0 * cont2csize + 4] = current_source[soffset + 112];
          current_target[toffset + 1 * cont2csize + 4] = current_source[soffset + 113];
          current_target[toffset + 2 * cont2csize + 4] = current_source[soffset + 114];
          current_target[toffset + 3 * cont2csize + 4] = current_source[soffset + 115];
          current_target[toffset + 4 * cont2csize + 4] = current_source[soffset + 116];
          current_target[toffset + 5 * cont2csize + 4] = current_source[soffset + 117];
          current_target[toffset + 6 * cont2csize + 4] = current_source[soffset + 118];
          current_target[toffset + 7 * cont2csize + 4] = current_source[soffset + 119];
          current_target[toffset + 8 * cont2csize + 4] = current_source[soffset + 120];
          current_target[toffset + 9 * cont2csize + 4] = current_source[soffset + 121];
          current_target[toffset + 10 * cont2csize + 4] = current_source[soffset + 122];
          current_target[toffset + 11 * cont2csize + 4] = current_source[soffset + 123];
          current_target[toffset + 12 * cont2csize + 4] = current_source[soffset + 124];
          current_target[toffset + 13 * cont2csize + 4] = current_source[soffset + 125];
          current_target[toffset + 14 * cont2csize + 4] = current_source[soffset + 126];
          current_target[toffset + 15 * cont2csize + 4] = current_source[soffset + 127];
          current_target[toffset + 16 * cont2csize + 4] = current_source[soffset + 128];
          current_target[toffset + 17 * cont2csize + 4] = current_source[soffset + 129];
          current_target[toffset + 18 * cont2csize + 4] = current_source[soffset + 130];
          current_target[toffset + 19 * cont2csize + 4] = current_source[soffset + 131];
          current_target[toffset + 20 * cont2csize + 4] = current_source[soffset + 132];
          current_target[toffset + 21 * cont2csize + 4] = current_source[soffset + 133];
          current_target[toffset + 22 * cont2csize + 4] = current_source[soffset + 134];
          current_target[toffset + 23 * cont2csize + 4] = current_source[soffset + 135];
          current_target[toffset + 24 * cont2csize + 4] = current_source[soffset + 136];
          current_target[toffset + 25 * cont2csize + 4] = current_source[soffset + 137];
          current_target[toffset + 26 * cont2csize + 4] = current_source[soffset + 138];
          current_target[toffset + 27 * cont2csize + 4] = current_source[soffset + 139];
          current_target[toffset + 0 * cont2csize + 5] = current_source[soffset + 140];
          current_target[toffset + 1 * cont2csize + 5] = current_source[soffset + 141];
          current_target[toffset + 2 * cont2csize + 5] = current_source[soffset + 142];
          current_target[toffset + 3 * cont2csize + 5] = current_source[soffset + 143];
          current_target[toffset + 4 * cont2csize + 5] = current_source[soffset + 144];
          current_target[toffset + 5 * cont2csize + 5] = current_source[soffset + 145];
          current_target[toffset + 6 * cont2csize + 5] = current_source[soffset + 146];
          current_target[toffset + 7 * cont2csize + 5] = current_source[soffset + 147];
          current_target[toffset + 8 * cont2csize + 5] = current_source[soffset + 148];
          current_target[toffset + 9 * cont2csize + 5] = current_source[soffset + 149];
          current_target[toffset + 10 * cont2csize + 5] = current_source[soffset + 150];
          current_target[toffset + 11 * cont2csize + 5] = current_source[soffset + 151];
          current_target[toffset + 12 * cont2csize + 5] = current_source[soffset + 152];
          current_target[toffset + 13 * cont2csize + 5] = current_source[soffset + 153];
          current_target[toffset + 14 * cont2csize + 5] = current_source[soffset + 154];
          current_target[toffset + 15 * cont2csize + 5] = current_source[soffset + 155];
          current_target[toffset + 16 * cont2csize + 5] = current_source[soffset + 156];
          current_target[toffset + 17 * cont2csize + 5] = current_source[soffset + 157];
          current_target[toffset + 18 * cont2csize + 5] = current_source[soffset + 158];
          current_target[toffset + 19 * cont2csize + 5] = current_source[soffset + 159];
          current_target[toffset + 20 * cont2csize + 5] = current_source[soffset + 160];
          current_target[toffset + 21 * cont2csize + 5] = current_source[soffset + 161];
          current_target[toffset + 22 * cont2csize + 5] = current_source[soffset + 162];
          current_target[toffset + 23 * cont2csize + 5] = current_source[soffset + 163];
          current_target[toffset + 24 * cont2csize + 5] = current_source[soffset + 164];
          current_target[toffset + 25 * cont2csize + 5] = current_source[soffset + 165];
          current_target[toffset + 26 * cont2csize + 5] = current_source[soffset + 166];
          current_target[toffset + 27 * cont2csize + 5] = current_source[soffset + 167];
          current_target[toffset + 0 * cont2csize + 6] = current_source[soffset + 168];
          current_target[toffset + 1 * cont2csize + 6] = current_source[soffset + 169];
          current_target[toffset + 2 * cont2csize + 6] = current_source[soffset + 170];
          current_target[toffset + 3 * cont2csize + 6] = current_source[soffset + 171];
          current_target[toffset + 4 * cont2csize + 6] = current_source[soffset + 172];
          current_target[toffset + 5 * cont2csize + 6] = current_source[soffset + 173];
          current_target[toffset + 6 * cont2csize + 6] = current_source[soffset + 174];
          current_target[toffset + 7 * cont2csize + 6] = current_source[soffset + 175];
          current_target[toffset + 8 * cont2csize + 6] = current_source[soffset + 176];
          current_target[toffset + 9 * cont2csize + 6] = current_source[soffset + 177];
          current_target[toffset + 10 * cont2csize + 6] = current_source[soffset + 178];
          current_target[toffset + 11 * cont2csize + 6] = current_source[soffset + 179];
          current_target[toffset + 12 * cont2csize + 6] = current_source[soffset + 180];
          current_target[toffset + 13 * cont2csize + 6] = current_source[soffset + 181];
          current_target[toffset + 14 * cont2csize + 6] = current_source[soffset + 182];
          current_target[toffset + 15 * cont2csize + 6] = current_source[soffset + 183];
          current_target[toffset + 16 * cont2csize + 6] = current_source[soffset + 184];
          current_target[toffset + 17 * cont2csize + 6] = current_source[soffset + 185];
          current_target[toffset + 18 * cont2csize + 6] = current_source[soffset + 186];
          current_target[toffset + 19 * cont2csize + 6] = current_source[soffset + 187];
          current_target[toffset + 20 * cont2csize + 6] = current_source[soffset + 188];
          current_target[toffset + 21 * cont2csize + 6] = current_source[soffset + 189];
          current_target[toffset + 22 * cont2csize + 6] = current_source[soffset + 190];
          current_target[toffset + 23 * cont2csize + 6] = current_source[soffset + 191];
          current_target[toffset + 24 * cont2csize + 6] = current_source[soffset + 192];
          current_target[toffset + 25 * cont2csize + 6] = current_source[soffset + 193];
          current_target[toffset + 26 * cont2csize + 6] = current_source[soffset + 194];
          current_target[toffset + 27 * cont2csize + 6] = current_source[soffset + 195];
          current_target[toffset + 0 * cont2csize + 7] = current_source[soffset + 196];
          current_target[toffset + 1 * cont2csize + 7] = current_source[soffset + 197];
          current_target[toffset + 2 * cont2csize + 7] = current_source[soffset + 198];
          current_target[toffset + 3 * cont2csize + 7] = current_source[soffset + 199];
          current_target[toffset + 4 * cont2csize + 7] = current_source[soffset + 200];
          current_target[toffset + 5 * cont2csize + 7] = current_source[soffset + 201];
          current_target[toffset + 6 * cont2csize + 7] = current_source[soffset + 202];
          current_target[toffset + 7 * cont2csize + 7] = current_source[soffset + 203];
          current_target[toffset + 8 * cont2csize + 7] = current_source[soffset + 204];
          current_target[toffset + 9 * cont2csize + 7] = current_source[soffset + 205];
          current_target[toffset + 10 * cont2csize + 7] = current_source[soffset + 206];
          current_target[toffset + 11 * cont2csize + 7] = current_source[soffset + 207];
          current_target[toffset + 12 * cont2csize + 7] = current_source[soffset + 208];
          current_target[toffset + 13 * cont2csize + 7] = current_source[soffset + 209];
          current_target[toffset + 14 * cont2csize + 7] = current_source[soffset + 210];
          current_target[toffset + 15 * cont2csize + 7] = current_source[soffset + 211];
          current_target[toffset + 16 * cont2csize + 7] = current_source[soffset + 212];
          current_target[toffset + 17 * cont2csize + 7] = current_source[soffset + 213];
          current_target[toffset + 18 * cont2csize + 7] = current_source[soffset + 214];
          current_target[toffset + 19 * cont2csize + 7] = current_source[soffset + 215];
          current_target[toffset + 20 * cont2csize + 7] = current_source[soffset + 216];
          current_target[toffset + 21 * cont2csize + 7] = current_source[soffset + 217];
          current_target[toffset + 22 * cont2csize + 7] = current_source[soffset + 218];
          current_target[toffset + 23 * cont2csize + 7] = current_source[soffset + 219];
          current_target[toffset + 24 * cont2csize + 7] = current_source[soffset + 220];
          current_target[toffset + 25 * cont2csize + 7] = current_source[soffset + 221];
          current_target[toffset + 26 * cont2csize + 7] = current_source[soffset + 222];
          current_target[toffset + 27 * cont2csize + 7] = current_source[soffset + 223];
          current_target[toffset + 0 * cont2csize + 8] = current_source[soffset + 224];
          current_target[toffset + 1 * cont2csize + 8] = current_source[soffset + 225];
          current_target[toffset + 2 * cont2csize + 8] = current_source[soffset + 226];
          current_target[toffset + 3 * cont2csize + 8] = current_source[soffset + 227];
          current_target[toffset + 4 * cont2csize + 8] = current_source[soffset + 228];
          current_target[toffset + 5 * cont2csize + 8] = current_source[soffset + 229];
          current_target[toffset + 6 * cont2csize + 8] = current_source[soffset + 230];
          current_target[toffset + 7 * cont2csize + 8] = current_source[soffset + 231];
          current_target[toffset + 8 * cont2csize + 8] = current_source[soffset + 232];
          current_target[toffset + 9 * cont2csize + 8] = current_source[soffset + 233];
          current_target[toffset + 10 * cont2csize + 8] = current_source[soffset + 234];
          current_target[toffset + 11 * cont2csize + 8] = current_source[soffset + 235];
          current_target[toffset + 12 * cont2csize + 8] = current_source[soffset + 236];
          current_target[toffset + 13 * cont2csize + 8] = current_source[soffset + 237];
          current_target[toffset + 14 * cont2csize + 8] = current_source[soffset + 238];
          current_target[toffset + 15 * cont2csize + 8] = current_source[soffset + 239];
          current_target[toffset + 16 * cont2csize + 8] = current_source[soffset + 240];
          current_target[toffset + 17 * cont2csize + 8] = current_source[soffset + 241];
          current_target[toffset + 18 * cont2csize + 8] = current_source[soffset + 242];
          current_target[toffset + 19 * cont2csize + 8] = current_source[soffset + 243];
          current_target[toffset + 20 * cont2csize + 8] = current_source[soffset + 244];
          current_target[toffset + 21 * cont2csize + 8] = current_source[soffset + 245];
          current_target[toffset + 22 * cont2csize + 8] = current_source[soffset + 246];
          current_target[toffset + 23 * cont2csize + 8] = current_source[soffset + 247];
          current_target[toffset + 24 * cont2csize + 8] = current_source[soffset + 248];
          current_target[toffset + 25 * cont2csize + 8] = current_source[soffset + 249];
          current_target[toffset + 26 * cont2csize + 8] = current_source[soffset + 250];
          current_target[toffset + 27 * cont2csize + 8] = current_source[soffset + 251];
          current_target[toffset + 0 * cont2csize + 9] = current_source[soffset + 252];
          current_target[toffset + 1 * cont2csize + 9] = current_source[soffset + 253];
          current_target[toffset + 2 * cont2csize + 9] = current_source[soffset + 254];
          current_target[toffset + 3 * cont2csize + 9] = current_source[soffset + 255];
          current_target[toffset + 4 * cont2csize + 9] = current_source[soffset + 256];
          current_target[toffset + 5 * cont2csize + 9] = current_source[soffset + 257];
          current_target[toffset + 6 * cont2csize + 9] = current_source[soffset + 258];
          current_target[toffset + 7 * cont2csize + 9] = current_source[soffset + 259];
          current_target[toffset + 8 * cont2csize + 9] = current_source[soffset + 260];
          current_target[toffset + 9 * cont2csize + 9] = current_source[soffset + 261];
          current_target[toffset + 10 * cont2csize + 9] = current_source[soffset + 262];
          current_target[toffset + 11 * cont2csize + 9] = current_source[soffset + 263];
          current_target[toffset + 12 * cont2csize + 9] = current_source[soffset + 264];
          current_target[toffset + 13 * cont2csize + 9] = current_source[soffset + 265];
          current_target[toffset + 14 * cont2csize + 9] = current_source[soffset + 266];
          current_target[toffset + 15 * cont2csize + 9] = current_source[soffset + 267];
          current_target[toffset + 16 * cont2csize + 9] = current_source[soffset + 268];
          current_target[toffset + 17 * cont2csize + 9] = current_source[soffset + 269];
          current_target[toffset + 18 * cont2csize + 9] = current_source[soffset + 270];
          current_target[toffset + 19 * cont2csize + 9] = current_source[soffset + 271];
          current_target[toffset + 20 * cont2csize + 9] = current_source[soffset + 272];
          current_target[toffset + 21 * cont2csize + 9] = current_source[soffset + 273];
          current_target[toffset + 22 * cont2csize + 9] = current_source[soffset + 274];
          current_target[toffset + 23 * cont2csize + 9] = current_source[soffset + 275];
          current_target[toffset + 24 * cont2csize + 9] = current_source[soffset + 276];
          current_target[toffset + 25 * cont2csize + 9] = current_source[soffset + 277];
          current_target[toffset + 26 * cont2csize + 9] = current_source[soffset + 278];
          current_target[toffset + 27 * cont2csize + 9] = current_source[soffset + 279];
          current_target[toffset + 0 * cont2csize + 10] = current_source[soffset + 280];
          current_target[toffset + 1 * cont2csize + 10] = current_source[soffset + 281];
          current_target[toffset + 2 * cont2csize + 10] = current_source[soffset + 282];
          current_target[toffset + 3 * cont2csize + 10] = current_source[soffset + 283];
          current_target[toffset + 4 * cont2csize + 10] = current_source[soffset + 284];
          current_target[toffset + 5 * cont2csize + 10] = current_source[soffset + 285];
          current_target[toffset + 6 * cont2csize + 10] = current_source[soffset + 286];
          current_target[toffset + 7 * cont2csize + 10] = current_source[soffset + 287];
          current_target[toffset + 8 * cont2csize + 10] = current_source[soffset + 288];
          current_target[toffset + 9 * cont2csize + 10] = current_source[soffset + 289];
          current_target[toffset + 10 * cont2csize + 10] = current_source[soffset + 290];
          current_target[toffset + 11 * cont2csize + 10] = current_source[soffset + 291];
          current_target[toffset + 12 * cont2csize + 10] = current_source[soffset + 292];
          current_target[toffset + 13 * cont2csize + 10] = current_source[soffset + 293];
          current_target[toffset + 14 * cont2csize + 10] = current_source[soffset + 294];
          current_target[toffset + 15 * cont2csize + 10] = current_source[soffset + 295];
          current_target[toffset + 16 * cont2csize + 10] = current_source[soffset + 296];
          current_target[toffset + 17 * cont2csize + 10] = current_source[soffset + 297];
          current_target[toffset + 18 * cont2csize + 10] = current_source[soffset + 298];
          current_target[toffset + 19 * cont2csize + 10] = current_source[soffset + 299];
          current_target[toffset + 20 * cont2csize + 10] = current_source[soffset + 300];
          current_target[toffset + 21 * cont2csize + 10] = current_source[soffset + 301];
          current_target[toffset + 22 * cont2csize + 10] = current_source[soffset + 302];
          current_target[toffset + 23 * cont2csize + 10] = current_source[soffset + 303];
          current_target[toffset + 24 * cont2csize + 10] = current_source[soffset + 304];
          current_target[toffset + 25 * cont2csize + 10] = current_source[soffset + 305];
          current_target[toffset + 26 * cont2csize + 10] = current_source[soffset + 306];
          current_target[toffset + 27 * cont2csize + 10] = current_source[soffset + 307];
          current_target[toffset + 0 * cont2csize + 11] = current_source[soffset + 308];
          current_target[toffset + 1 * cont2csize + 11] = current_source[soffset + 309];
          current_target[toffset + 2 * cont2csize + 11] = current_source[soffset + 310];
          current_target[toffset + 3 * cont2csize + 11] = current_source[soffset + 311];
          current_target[toffset + 4 * cont2csize + 11] = current_source[soffset + 312];
          current_target[toffset + 5 * cont2csize + 11] = current_source[soffset + 313];
          current_target[toffset + 6 * cont2csize + 11] = current_source[soffset + 314];
          current_target[toffset + 7 * cont2csize + 11] = current_source[soffset + 315];
          current_target[toffset + 8 * cont2csize + 11] = current_source[soffset + 316];
          current_target[toffset + 9 * cont2csize + 11] = current_source[soffset + 317];
          current_target[toffset + 10 * cont2csize + 11] = current_source[soffset + 318];
          current_target[toffset + 11 * cont2csize + 11] = current_source[soffset + 319];
          current_target[toffset + 12 * cont2csize + 11] = current_source[soffset + 320];
          current_target[toffset + 13 * cont2csize + 11] = current_source[soffset + 321];
          current_target[toffset + 14 * cont2csize + 11] = current_source[soffset + 322];
          current_target[toffset + 15 * cont2csize + 11] = current_source[soffset + 323];
          current_target[toffset + 16 * cont2csize + 11] = current_source[soffset + 324];
          current_target[toffset + 17 * cont2csize + 11] = current_source[soffset + 325];
          current_target[toffset + 18 * cont2csize + 11] = current_source[soffset + 326];
          current_target[toffset + 19 * cont2csize + 11] = current_source[soffset + 327];
          current_target[toffset + 20 * cont2csize + 11] = current_source[soffset + 328];
          current_target[toffset + 21 * cont2csize + 11] = current_source[soffset + 329];
          current_target[toffset + 22 * cont2csize + 11] = current_source[soffset + 330];
          current_target[toffset + 23 * cont2csize + 11] = current_source[soffset + 331];
          current_target[toffset + 24 * cont2csize + 11] = current_source[soffset + 332];
          current_target[toffset + 25 * cont2csize + 11] = current_source[soffset + 333];
          current_target[toffset + 26 * cont2csize + 11] = current_source[soffset + 334];
          current_target[toffset + 27 * cont2csize + 11] = current_source[soffset + 335];
          current_target[toffset + 0 * cont2csize + 12] = current_source[soffset + 336];
          current_target[toffset + 1 * cont2csize + 12] = current_source[soffset + 337];
          current_target[toffset + 2 * cont2csize + 12] = current_source[soffset + 338];
          current_target[toffset + 3 * cont2csize + 12] = current_source[soffset + 339];
          current_target[toffset + 4 * cont2csize + 12] = current_source[soffset + 340];
          current_target[toffset + 5 * cont2csize + 12] = current_source[soffset + 341];
          current_target[toffset + 6 * cont2csize + 12] = current_source[soffset + 342];
          current_target[toffset + 7 * cont2csize + 12] = current_source[soffset + 343];
          current_target[toffset + 8 * cont2csize + 12] = current_source[soffset + 344];
          current_target[toffset + 9 * cont2csize + 12] = current_source[soffset + 345];
          current_target[toffset + 10 * cont2csize + 12] = current_source[soffset + 346];
          current_target[toffset + 11 * cont2csize + 12] = current_source[soffset + 347];
          current_target[toffset + 12 * cont2csize + 12] = current_source[soffset + 348];
          current_target[toffset + 13 * cont2csize + 12] = current_source[soffset + 349];
          current_target[toffset + 14 * cont2csize + 12] = current_source[soffset + 350];
          current_target[toffset + 15 * cont2csize + 12] = current_source[soffset + 351];
          current_target[toffset + 16 * cont2csize + 12] = current_source[soffset + 352];
          current_target[toffset + 17 * cont2csize + 12] = current_source[soffset + 353];
          current_target[toffset + 18 * cont2csize + 12] = current_source[soffset + 354];
          current_target[toffset + 19 * cont2csize + 12] = current_source[soffset + 355];
          current_target[toffset + 20 * cont2csize + 12] = current_source[soffset + 356];
          current_target[toffset + 21 * cont2csize + 12] = current_source[soffset + 357];
          current_target[toffset + 22 * cont2csize + 12] = current_source[soffset + 358];
          current_target[toffset + 23 * cont2csize + 12] = current_source[soffset + 359];
          current_target[toffset + 24 * cont2csize + 12] = current_source[soffset + 360];
          current_target[toffset + 25 * cont2csize + 12] = current_source[soffset + 361];
          current_target[toffset + 26 * cont2csize + 12] = current_source[soffset + 362];
          current_target[toffset + 27 * cont2csize + 12] = current_source[soffset + 363];
          current_target[toffset + 0 * cont2csize + 13] = current_source[soffset + 364];
          current_target[toffset + 1 * cont2csize + 13] = current_source[soffset + 365];
          current_target[toffset + 2 * cont2csize + 13] = current_source[soffset + 366];
          current_target[toffset + 3 * cont2csize + 13] = current_source[soffset + 367];
          current_target[toffset + 4 * cont2csize + 13] = current_source[soffset + 368];
          current_target[toffset + 5 * cont2csize + 13] = current_source[soffset + 369];
          current_target[toffset + 6 * cont2csize + 13] = current_source[soffset + 370];
          current_target[toffset + 7 * cont2csize + 13] = current_source[soffset + 371];
          current_target[toffset + 8 * cont2csize + 13] = current_source[soffset + 372];
          current_target[toffset + 9 * cont2csize + 13] = current_source[soffset + 373];
          current_target[toffset + 10 * cont2csize + 13] = current_source[soffset + 374];
          current_target[toffset + 11 * cont2csize + 13] = current_source[soffset + 375];
          current_target[toffset + 12 * cont2csize + 13] = current_source[soffset + 376];
          current_target[toffset + 13 * cont2csize + 13] = current_source[soffset + 377];
          current_target[toffset + 14 * cont2csize + 13] = current_source[soffset + 378];
          current_target[toffset + 15 * cont2csize + 13] = current_source[soffset + 379];
          current_target[toffset + 16 * cont2csize + 13] = current_source[soffset + 380];
          current_target[toffset + 17 * cont2csize + 13] = current_source[soffset + 381];
          current_target[toffset + 18 * cont2csize + 13] = current_source[soffset + 382];
          current_target[toffset + 19 * cont2csize + 13] = current_source[soffset + 383];
          current_target[toffset + 20 * cont2csize + 13] = current_source[soffset + 384];
          current_target[toffset + 21 * cont2csize + 13] = current_source[soffset + 385];
          current_target[toffset + 22 * cont2csize + 13] = current_source[soffset + 386];
          current_target[toffset + 23 * cont2csize + 13] = current_source[soffset + 387];
          current_target[toffset + 24 * cont2csize + 13] = current_source[soffset + 388];
          current_target[toffset + 25 * cont2csize + 13] = current_source[soffset + 389];
          current_target[toffset + 26 * cont2csize + 13] = current_source[soffset + 390];
          current_target[toffset + 27 * cont2csize + 13] = current_source[soffset + 391];
          current_target[toffset + 0 * cont2csize + 14] = current_source[soffset + 392];
          current_target[toffset + 1 * cont2csize + 14] = current_source[soffset + 393];
          current_target[toffset + 2 * cont2csize + 14] = current_source[soffset + 394];
          current_target[toffset + 3 * cont2csize + 14] = current_source[soffset + 395];
          current_target[toffset + 4 * cont2csize + 14] = current_source[soffset + 396];
          current_target[toffset + 5 * cont2csize + 14] = current_source[soffset + 397];
          current_target[toffset + 6 * cont2csize + 14] = current_source[soffset + 398];
          current_target[toffset + 7 * cont2csize + 14] = current_source[soffset + 399];
          current_target[toffset + 8 * cont2csize + 14] = current_source[soffset + 400];
          current_target[toffset + 9 * cont2csize + 14] = current_source[soffset + 401];
          current_target[toffset + 10 * cont2csize + 14] = current_source[soffset + 402];
          current_target[toffset + 11 * cont2csize + 14] = current_source[soffset + 403];
          current_target[toffset + 12 * cont2csize + 14] = current_source[soffset + 404];
          current_target[toffset + 13 * cont2csize + 14] = current_source[soffset + 405];
          current_target[toffset + 14 * cont2csize + 14] = current_source[soffset + 406];
          current_target[toffset + 15 * cont2csize + 14] = current_source[soffset + 407];
          current_target[toffset + 16 * cont2csize + 14] = current_source[soffset + 408];
          current_target[toffset + 17 * cont2csize + 14] = current_source[soffset + 409];
          current_target[toffset + 18 * cont2csize + 14] = current_source[soffset + 410];
          current_target[toffset + 19 * cont2csize + 14] = current_source[soffset + 411];
          current_target[toffset + 20 * cont2csize + 14] = current_source[soffset + 412];
          current_target[toffset + 21 * cont2csize + 14] = current_source[soffset + 413];
          current_target[toffset + 22 * cont2csize + 14] = current_source[soffset + 414];
          current_target[toffset + 23 * cont2csize + 14] = current_source[soffset + 415];
          current_target[toffset + 24 * cont2csize + 14] = current_source[soffset + 416];
          current_target[toffset + 25 * cont2csize + 14] = current_source[soffset + 417];
          current_target[toffset + 26 * cont2csize + 14] = current_source[soffset + 418];
          current_target[toffset + 27 * cont2csize + 14] = current_source[soffset + 419];
          current_target[toffset + 0 * cont2csize + 15] = current_source[soffset + 420];
          current_target[toffset + 1 * cont2csize + 15] = current_source[soffset + 421];
          current_target[toffset + 2 * cont2csize + 15] = current_source[soffset + 422];
          current_target[toffset + 3 * cont2csize + 15] = current_source[soffset + 423];
          current_target[toffset + 4 * cont2csize + 15] = current_source[soffset + 424];
          current_target[toffset + 5 * cont2csize + 15] = current_source[soffset + 425];
          current_target[toffset + 6 * cont2csize + 15] = current_source[soffset + 426];
          current_target[toffset + 7 * cont2csize + 15] = current_source[soffset + 427];
          current_target[toffset + 8 * cont2csize + 15] = current_source[soffset + 428];
          current_target[toffset + 9 * cont2csize + 15] = current_source[soffset + 429];
          current_target[toffset + 10 * cont2csize + 15] = current_source[soffset + 430];
          current_target[toffset + 11 * cont2csize + 15] = current_source[soffset + 431];
          current_target[toffset + 12 * cont2csize + 15] = current_source[soffset + 432];
          current_target[toffset + 13 * cont2csize + 15] = current_source[soffset + 433];
          current_target[toffset + 14 * cont2csize + 15] = current_source[soffset + 434];
          current_target[toffset + 15 * cont2csize + 15] = current_source[soffset + 435];
          current_target[toffset + 16 * cont2csize + 15] = current_source[soffset + 436];
          current_target[toffset + 17 * cont2csize + 15] = current_source[soffset + 437];
          current_target[toffset + 18 * cont2csize + 15] = current_source[soffset + 438];
          current_target[toffset + 19 * cont2csize + 15] = current_source[soffset + 439];
          current_target[toffset + 20 * cont2csize + 15] = current_source[soffset + 440];
          current_target[toffset + 21 * cont2csize + 15] = current_source[soffset + 441];
          current_target[toffset + 22 * cont2csize + 15] = current_source[soffset + 442];
          current_target[toffset + 23 * cont2csize + 15] = current_source[soffset + 443];
          current_target[toffset + 24 * cont2csize + 15] = current_source[soffset + 444];
          current_target[toffset + 25 * cont2csize + 15] = current_source[soffset + 445];
          current_target[toffset + 26 * cont2csize + 15] = current_source[soffset + 446];
          current_target[toffset + 27 * cont2csize + 15] = current_source[soffset + 447];
          current_target[toffset + 0 * cont2csize + 16] = current_source[soffset + 448];
          current_target[toffset + 1 * cont2csize + 16] = current_source[soffset + 449];
          current_target[toffset + 2 * cont2csize + 16] = current_source[soffset + 450];
          current_target[toffset + 3 * cont2csize + 16] = current_source[soffset + 451];
          current_target[toffset + 4 * cont2csize + 16] = current_source[soffset + 452];
          current_target[toffset + 5 * cont2csize + 16] = current_source[soffset + 453];
          current_target[toffset + 6 * cont2csize + 16] = current_source[soffset + 454];
          current_target[toffset + 7 * cont2csize + 16] = current_source[soffset + 455];
          current_target[toffset + 8 * cont2csize + 16] = current_source[soffset + 456];
          current_target[toffset + 9 * cont2csize + 16] = current_source[soffset + 457];
          current_target[toffset + 10 * cont2csize + 16] = current_source[soffset + 458];
          current_target[toffset + 11 * cont2csize + 16] = current_source[soffset + 459];
          current_target[toffset + 12 * cont2csize + 16] = current_source[soffset + 460];
          current_target[toffset + 13 * cont2csize + 16] = current_source[soffset + 461];
          current_target[toffset + 14 * cont2csize + 16] = current_source[soffset + 462];
          current_target[toffset + 15 * cont2csize + 16] = current_source[soffset + 463];
          current_target[toffset + 16 * cont2csize + 16] = current_source[soffset + 464];
          current_target[toffset + 17 * cont2csize + 16] = current_source[soffset + 465];
          current_target[toffset + 18 * cont2csize + 16] = current_source[soffset + 466];
          current_target[toffset + 19 * cont2csize + 16] = current_source[soffset + 467];
          current_target[toffset + 20 * cont2csize + 16] = current_source[soffset + 468];
          current_target[toffset + 21 * cont2csize + 16] = current_source[soffset + 469];
          current_target[toffset + 22 * cont2csize + 16] = current_source[soffset + 470];
          current_target[toffset + 23 * cont2csize + 16] = current_source[soffset + 471];
          current_target[toffset + 24 * cont2csize + 16] = current_source[soffset + 472];
          current_target[toffset + 25 * cont2csize + 16] = current_source[soffset + 473];
          current_target[toffset + 26 * cont2csize + 16] = current_source[soffset + 474];
          current_target[toffset + 27 * cont2csize + 16] = current_source[soffset + 475];
          current_target[toffset + 0 * cont2csize + 17] = current_source[soffset + 476];
          current_target[toffset + 1 * cont2csize + 17] = current_source[soffset + 477];
          current_target[toffset + 2 * cont2csize + 17] = current_source[soffset + 478];
          current_target[toffset + 3 * cont2csize + 17] = current_source[soffset + 479];
          current_target[toffset + 4 * cont2csize + 17] = current_source[soffset + 480];
          current_target[toffset + 5 * cont2csize + 17] = current_source[soffset + 481];
          current_target[toffset + 6 * cont2csize + 17] = current_source[soffset + 482];
          current_target[toffset + 7 * cont2csize + 17] = current_source[soffset + 483];
          current_target[toffset + 8 * cont2csize + 17] = current_source[soffset + 484];
          current_target[toffset + 9 * cont2csize + 17] = current_source[soffset + 485];
          current_target[toffset + 10 * cont2csize + 17] = current_source[soffset + 486];
          current_target[toffset + 11 * cont2csize + 17] = current_source[soffset + 487];
          current_target[toffset + 12 * cont2csize + 17] = current_source[soffset + 488];
          current_target[toffset + 13 * cont2csize + 17] = current_source[soffset + 489];
          current_target[toffset + 14 * cont2csize + 17] = current_source[soffset + 490];
          current_target[toffset + 15 * cont2csize + 17] = current_source[soffset + 491];
          current_target[toffset + 16 * cont2csize + 17] = current_source[soffset + 492];
          current_target[toffset + 17 * cont2csize + 17] = current_source[soffset + 493];
          current_target[toffset + 18 * cont2csize + 17] = current_source[soffset + 494];
          current_target[toffset + 19 * cont2csize + 17] = current_source[soffset + 495];
          current_target[toffset + 20 * cont2csize + 17] = current_source[soffset + 496];
          current_target[toffset + 21 * cont2csize + 17] = current_source[soffset + 497];
          current_target[toffset + 22 * cont2csize + 17] = current_source[soffset + 498];
          current_target[toffset + 23 * cont2csize + 17] = current_source[soffset + 499];
          current_target[toffset + 24 * cont2csize + 17] = current_source[soffset + 500];
          current_target[toffset + 25 * cont2csize + 17] = current_source[soffset + 501];
          current_target[toffset + 26 * cont2csize + 17] = current_source[soffset + 502];
          current_target[toffset + 27 * cont2csize + 17] = current_source[soffset + 503];
          current_target[toffset + 0 * cont2csize + 18] = current_source[soffset + 504];
          current_target[toffset + 1 * cont2csize + 18] = current_source[soffset + 505];
          current_target[toffset + 2 * cont2csize + 18] = current_source[soffset + 506];
          current_target[toffset + 3 * cont2csize + 18] = current_source[soffset + 507];
          current_target[toffset + 4 * cont2csize + 18] = current_source[soffset + 508];
          current_target[toffset + 5 * cont2csize + 18] = current_source[soffset + 509];
          current_target[toffset + 6 * cont2csize + 18] = current_source[soffset + 510];
          current_target[toffset + 7 * cont2csize + 18] = current_source[soffset + 511];
          current_target[toffset + 8 * cont2csize + 18] = current_source[soffset + 512];
          current_target[toffset + 9 * cont2csize + 18] = current_source[soffset + 513];
          current_target[toffset + 10 * cont2csize + 18] = current_source[soffset + 514];
          current_target[toffset + 11 * cont2csize + 18] = current_source[soffset + 515];
          current_target[toffset + 12 * cont2csize + 18] = current_source[soffset + 516];
          current_target[toffset + 13 * cont2csize + 18] = current_source[soffset + 517];
          current_target[toffset + 14 * cont2csize + 18] = current_source[soffset + 518];
          current_target[toffset + 15 * cont2csize + 18] = current_source[soffset + 519];
          current_target[toffset + 16 * cont2csize + 18] = current_source[soffset + 520];
          current_target[toffset + 17 * cont2csize + 18] = current_source[soffset + 521];
          current_target[toffset + 18 * cont2csize + 18] = current_source[soffset + 522];
          current_target[toffset + 19 * cont2csize + 18] = current_source[soffset + 523];
          current_target[toffset + 20 * cont2csize + 18] = current_source[soffset + 524];
          current_target[toffset + 21 * cont2csize + 18] = current_source[soffset + 525];
          current_target[toffset + 22 * cont2csize + 18] = current_source[soffset + 526];
          current_target[toffset + 23 * cont2csize + 18] = current_source[soffset + 527];
          current_target[toffset + 24 * cont2csize + 18] = current_source[soffset + 528];
          current_target[toffset + 25 * cont2csize + 18] = current_source[soffset + 529];
          current_target[toffset + 26 * cont2csize + 18] = current_source[soffset + 530];
          current_target[toffset + 27 * cont2csize + 18] = current_source[soffset + 531];
          current_target[toffset + 0 * cont2csize + 19] = current_source[soffset + 532];
          current_target[toffset + 1 * cont2csize + 19] = current_source[soffset + 533];
          current_target[toffset + 2 * cont2csize + 19] = current_source[soffset + 534];
          current_target[toffset + 3 * cont2csize + 19] = current_source[soffset + 535];
          current_target[toffset + 4 * cont2csize + 19] = current_source[soffset + 536];
          current_target[toffset + 5 * cont2csize + 19] = current_source[soffset + 537];
          current_target[toffset + 6 * cont2csize + 19] = current_source[soffset + 538];
          current_target[toffset + 7 * cont2csize + 19] = current_source[soffset + 539];
          current_target[toffset + 8 * cont2csize + 19] = current_source[soffset + 540];
          current_target[toffset + 9 * cont2csize + 19] = current_source[soffset + 541];
          current_target[toffset + 10 * cont2csize + 19] = current_source[soffset + 542];
          current_target[toffset + 11 * cont2csize + 19] = current_source[soffset + 543];
          current_target[toffset + 12 * cont2csize + 19] = current_source[soffset + 544];
          current_target[toffset + 13 * cont2csize + 19] = current_source[soffset + 545];
          current_target[toffset + 14 * cont2csize + 19] = current_source[soffset + 546];
          current_target[toffset + 15 * cont2csize + 19] = current_source[soffset + 547];
          current_target[toffset + 16 * cont2csize + 19] = current_source[soffset + 548];
          current_target[toffset + 17 * cont2csize + 19] = current_source[soffset + 549];
          current_target[toffset + 18 * cont2csize + 19] = current_source[soffset + 550];
          current_target[toffset + 19 * cont2csize + 19] = current_source[soffset + 551];
          current_target[toffset + 20 * cont2csize + 19] = current_source[soffset + 552];
          current_target[toffset + 21 * cont2csize + 19] = current_source[soffset + 553];
          current_target[toffset + 22 * cont2csize + 19] = current_source[soffset + 554];
          current_target[toffset + 23 * cont2csize + 19] = current_source[soffset + 555];
          current_target[toffset + 24 * cont2csize + 19] = current_source[soffset + 556];
          current_target[toffset + 25 * cont2csize + 19] = current_source[soffset + 557];
          current_target[toffset + 26 * cont2csize + 19] = current_source[soffset + 558];
          current_target[toffset + 27 * cont2csize + 19] = current_source[soffset + 559];
          current_target[toffset + 0 * cont2csize + 20] = current_source[soffset + 560];
          current_target[toffset + 1 * cont2csize + 20] = current_source[soffset + 561];
          current_target[toffset + 2 * cont2csize + 20] = current_source[soffset + 562];
          current_target[toffset + 3 * cont2csize + 20] = current_source[soffset + 563];
          current_target[toffset + 4 * cont2csize + 20] = current_source[soffset + 564];
          current_target[toffset + 5 * cont2csize + 20] = current_source[soffset + 565];
          current_target[toffset + 6 * cont2csize + 20] = current_source[soffset + 566];
          current_target[toffset + 7 * cont2csize + 20] = current_source[soffset + 567];
          current_target[toffset + 8 * cont2csize + 20] = current_source[soffset + 568];
          current_target[toffset + 9 * cont2csize + 20] = current_source[soffset + 569];
          current_target[toffset + 10 * cont2csize + 20] = current_source[soffset + 570];
          current_target[toffset + 11 * cont2csize + 20] = current_source[soffset + 571];
          current_target[toffset + 12 * cont2csize + 20] = current_source[soffset + 572];
          current_target[toffset + 13 * cont2csize + 20] = current_source[soffset + 573];
          current_target[toffset + 14 * cont2csize + 20] = current_source[soffset + 574];
          current_target[toffset + 15 * cont2csize + 20] = current_source[soffset + 575];
          current_target[toffset + 16 * cont2csize + 20] = current_source[soffset + 576];
          current_target[toffset + 17 * cont2csize + 20] = current_source[soffset + 577];
          current_target[toffset + 18 * cont2csize + 20] = current_source[soffset + 578];
          current_target[toffset + 19 * cont2csize + 20] = current_source[soffset + 579];
          current_target[toffset + 20 * cont2csize + 20] = current_source[soffset + 580];
          current_target[toffset + 21 * cont2csize + 20] = current_source[soffset + 581];
          current_target[toffset + 22 * cont2csize + 20] = current_source[soffset + 582];
          current_target[toffset + 23 * cont2csize + 20] = current_source[soffset + 583];
          current_target[toffset + 24 * cont2csize + 20] = current_source[soffset + 584];
          current_target[toffset + 25 * cont2csize + 20] = current_source[soffset + 585];
          current_target[toffset + 26 * cont2csize + 20] = current_source[soffset + 586];
          current_target[toffset + 27 * cont2csize + 20] = current_source[soffset + 587];
          current_target[toffset + 0 * cont2csize + 21] = current_source[soffset + 588];
          current_target[toffset + 1 * cont2csize + 21] = current_source[soffset + 589];
          current_target[toffset + 2 * cont2csize + 21] = current_source[soffset + 590];
          current_target[toffset + 3 * cont2csize + 21] = current_source[soffset + 591];
          current_target[toffset + 4 * cont2csize + 21] = current_source[soffset + 592];
          current_target[toffset + 5 * cont2csize + 21] = current_source[soffset + 593];
          current_target[toffset + 6 * cont2csize + 21] = current_source[soffset + 594];
          current_target[toffset + 7 * cont2csize + 21] = current_source[soffset + 595];
          current_target[toffset + 8 * cont2csize + 21] = current_source[soffset + 596];
          current_target[toffset + 9 * cont2csize + 21] = current_source[soffset + 597];
          current_target[toffset + 10 * cont2csize + 21] = current_source[soffset + 598];
          current_target[toffset + 11 * cont2csize + 21] = current_source[soffset + 599];
          current_target[toffset + 12 * cont2csize + 21] = current_source[soffset + 600];
          current_target[toffset + 13 * cont2csize + 21] = current_source[soffset + 601];
          current_target[toffset + 14 * cont2csize + 21] = current_source[soffset + 602];
          current_target[toffset + 15 * cont2csize + 21] = current_source[soffset + 603];
          current_target[toffset + 16 * cont2csize + 21] = current_source[soffset + 604];
          current_target[toffset + 17 * cont2csize + 21] = current_source[soffset + 605];
          current_target[toffset + 18 * cont2csize + 21] = current_source[soffset + 606];
          current_target[toffset + 19 * cont2csize + 21] = current_source[soffset + 607];
          current_target[toffset + 20 * cont2csize + 21] = current_source[soffset + 608];
          current_target[toffset + 21 * cont2csize + 21] = current_source[soffset + 609];
          current_target[toffset + 22 * cont2csize + 21] = current_source[soffset + 610];
          current_target[toffset + 23 * cont2csize + 21] = current_source[soffset + 611];
          current_target[toffset + 24 * cont2csize + 21] = current_source[soffset + 612];
          current_target[toffset + 25 * cont2csize + 21] = current_source[soffset + 613];
          current_target[toffset + 26 * cont2csize + 21] = current_source[soffset + 614];
          current_target[toffset + 27 * cont2csize + 21] = current_source[soffset + 615];
          current_target[toffset + 0 * cont2csize + 22] = current_source[soffset + 616];
          current_target[toffset + 1 * cont2csize + 22] = current_source[soffset + 617];
          current_target[toffset + 2 * cont2csize + 22] = current_source[soffset + 618];
          current_target[toffset + 3 * cont2csize + 22] = current_source[soffset + 619];
          current_target[toffset + 4 * cont2csize + 22] = current_source[soffset + 620];
          current_target[toffset + 5 * cont2csize + 22] = current_source[soffset + 621];
          current_target[toffset + 6 * cont2csize + 22] = current_source[soffset + 622];
          current_target[toffset + 7 * cont2csize + 22] = current_source[soffset + 623];
          current_target[toffset + 8 * cont2csize + 22] = current_source[soffset + 624];
          current_target[toffset + 9 * cont2csize + 22] = current_source[soffset + 625];
          current_target[toffset + 10 * cont2csize + 22] = current_source[soffset + 626];
          current_target[toffset + 11 * cont2csize + 22] = current_source[soffset + 627];
          current_target[toffset + 12 * cont2csize + 22] = current_source[soffset + 628];
          current_target[toffset + 13 * cont2csize + 22] = current_source[soffset + 629];
          current_target[toffset + 14 * cont2csize + 22] = current_source[soffset + 630];
          current_target[toffset + 15 * cont2csize + 22] = current_source[soffset + 631];
          current_target[toffset + 16 * cont2csize + 22] = current_source[soffset + 632];
          current_target[toffset + 17 * cont2csize + 22] = current_source[soffset + 633];
          current_target[toffset + 18 * cont2csize + 22] = current_source[soffset + 634];
          current_target[toffset + 19 * cont2csize + 22] = current_source[soffset + 635];
          current_target[toffset + 20 * cont2csize + 22] = current_source[soffset + 636];
          current_target[toffset + 21 * cont2csize + 22] = current_source[soffset + 637];
          current_target[toffset + 22 * cont2csize + 22] = current_source[soffset + 638];
          current_target[toffset + 23 * cont2csize + 22] = current_source[soffset + 639];
          current_target[toffset + 24 * cont2csize + 22] = current_source[soffset + 640];
          current_target[toffset + 25 * cont2csize + 22] = current_source[soffset + 641];
          current_target[toffset + 26 * cont2csize + 22] = current_source[soffset + 642];
          current_target[toffset + 27 * cont2csize + 22] = current_source[soffset + 643];
          current_target[toffset + 0 * cont2csize + 23] = current_source[soffset + 644];
          current_target[toffset + 1 * cont2csize + 23] = current_source[soffset + 645];
          current_target[toffset + 2 * cont2csize + 23] = current_source[soffset + 646];
          current_target[toffset + 3 * cont2csize + 23] = current_source[soffset + 647];
          current_target[toffset + 4 * cont2csize + 23] = current_source[soffset + 648];
          current_target[toffset + 5 * cont2csize + 23] = current_source[soffset + 649];
          current_target[toffset + 6 * cont2csize + 23] = current_source[soffset + 650];
          current_target[toffset + 7 * cont2csize + 23] = current_source[soffset + 651];
          current_target[toffset + 8 * cont2csize + 23] = current_source[soffset + 652];
          current_target[toffset + 9 * cont2csize + 23] = current_source[soffset + 653];
          current_target[toffset + 10 * cont2csize + 23] = current_source[soffset + 654];
          current_target[toffset + 11 * cont2csize + 23] = current_source[soffset + 655];
          current_target[toffset + 12 * cont2csize + 23] = current_source[soffset + 656];
          current_target[toffset + 13 * cont2csize + 23] = current_source[soffset + 657];
          current_target[toffset + 14 * cont2csize + 23] = current_source[soffset + 658];
          current_target[toffset + 15 * cont2csize + 23] = current_source[soffset + 659];
          current_target[toffset + 16 * cont2csize + 23] = current_source[soffset + 660];
          current_target[toffset + 17 * cont2csize + 23] = current_source[soffset + 661];
          current_target[toffset + 18 * cont2csize + 23] = current_source[soffset + 662];
          current_target[toffset + 19 * cont2csize + 23] = current_source[soffset + 663];
          current_target[toffset + 20 * cont2csize + 23] = current_source[soffset + 664];
          current_target[toffset + 21 * cont2csize + 23] = current_source[soffset + 665];
          current_target[toffset + 22 * cont2csize + 23] = current_source[soffset + 666];
          current_target[toffset + 23 * cont2csize + 23] = current_source[soffset + 667];
          current_target[toffset + 24 * cont2csize + 23] = current_source[soffset + 668];
          current_target[toffset + 25 * cont2csize + 23] = current_source[soffset + 669];
          current_target[toffset + 26 * cont2csize + 23] = current_source[soffset + 670];
          current_target[toffset + 27 * cont2csize + 23] = current_source[soffset + 671];
          current_target[toffset + 0 * cont2csize + 24] = current_source[soffset + 672];
          current_target[toffset + 1 * cont2csize + 24] = current_source[soffset + 673];
          current_target[toffset + 2 * cont2csize + 24] = current_source[soffset + 674];
          current_target[toffset + 3 * cont2csize + 24] = current_source[soffset + 675];
          current_target[toffset + 4 * cont2csize + 24] = current_source[soffset + 676];
          current_target[toffset + 5 * cont2csize + 24] = current_source[soffset + 677];
          current_target[toffset + 6 * cont2csize + 24] = current_source[soffset + 678];
          current_target[toffset + 7 * cont2csize + 24] = current_source[soffset + 679];
          current_target[toffset + 8 * cont2csize + 24] = current_source[soffset + 680];
          current_target[toffset + 9 * cont2csize + 24] = current_source[soffset + 681];
          current_target[toffset + 10 * cont2csize + 24] = current_source[soffset + 682];
          current_target[toffset + 11 * cont2csize + 24] = current_source[soffset + 683];
          current_target[toffset + 12 * cont2csize + 24] = current_source[soffset + 684];
          current_target[toffset + 13 * cont2csize + 24] = current_source[soffset + 685];
          current_target[toffset + 14 * cont2csize + 24] = current_source[soffset + 686];
          current_target[toffset + 15 * cont2csize + 24] = current_source[soffset + 687];
          current_target[toffset + 16 * cont2csize + 24] = current_source[soffset + 688];
          current_target[toffset + 17 * cont2csize + 24] = current_source[soffset + 689];
          current_target[toffset + 18 * cont2csize + 24] = current_source[soffset + 690];
          current_target[toffset + 19 * cont2csize + 24] = current_source[soffset + 691];
          current_target[toffset + 20 * cont2csize + 24] = current_source[soffset + 692];
          current_target[toffset + 21 * cont2csize + 24] = current_source[soffset + 693];
          current_target[toffset + 22 * cont2csize + 24] = current_source[soffset + 694];
          current_target[toffset + 23 * cont2csize + 24] = current_source[soffset + 695];
          current_target[toffset + 24 * cont2csize + 24] = current_source[soffset + 696];
          current_target[toffset + 25 * cont2csize + 24] = current_source[soffset + 697];
          current_target[toffset + 26 * cont2csize + 24] = current_source[soffset + 698];
          current_target[toffset + 27 * cont2csize + 24] = current_source[soffset + 699];
          current_target[toffset + 0 * cont2csize + 25] = current_source[soffset + 700];
          current_target[toffset + 1 * cont2csize + 25] = current_source[soffset + 701];
          current_target[toffset + 2 * cont2csize + 25] = current_source[soffset + 702];
          current_target[toffset + 3 * cont2csize + 25] = current_source[soffset + 703];
          current_target[toffset + 4 * cont2csize + 25] = current_source[soffset + 704];
          current_target[toffset + 5 * cont2csize + 25] = current_source[soffset + 705];
          current_target[toffset + 6 * cont2csize + 25] = current_source[soffset + 706];
          current_target[toffset + 7 * cont2csize + 25] = current_source[soffset + 707];
          current_target[toffset + 8 * cont2csize + 25] = current_source[soffset + 708];
          current_target[toffset + 9 * cont2csize + 25] = current_source[soffset + 709];
          current_target[toffset + 10 * cont2csize + 25] = current_source[soffset + 710];
          current_target[toffset + 11 * cont2csize + 25] = current_source[soffset + 711];
          current_target[toffset + 12 * cont2csize + 25] = current_source[soffset + 712];
          current_target[toffset + 13 * cont2csize + 25] = current_source[soffset + 713];
          current_target[toffset + 14 * cont2csize + 25] = current_source[soffset + 714];
          current_target[toffset + 15 * cont2csize + 25] = current_source[soffset + 715];
          current_target[toffset + 16 * cont2csize + 25] = current_source[soffset + 716];
          current_target[toffset + 17 * cont2csize + 25] = current_source[soffset + 717];
          current_target[toffset + 18 * cont2csize + 25] = current_source[soffset + 718];
          current_target[toffset + 19 * cont2csize + 25] = current_source[soffset + 719];
          current_target[toffset + 20 * cont2csize + 25] = current_source[soffset + 720];
          current_target[toffset + 21 * cont2csize + 25] = current_source[soffset + 721];
          current_target[toffset + 22 * cont2csize + 25] = current_source[soffset + 722];
          current_target[toffset + 23 * cont2csize + 25] = current_source[soffset + 723];
          current_target[toffset + 24 * cont2csize + 25] = current_source[soffset + 724];
          current_target[toffset + 25 * cont2csize + 25] = current_source[soffset + 725];
          current_target[toffset + 26 * cont2csize + 25] = current_source[soffset + 726];
          current_target[toffset + 27 * cont2csize + 25] = current_source[soffset + 727];
          current_target[toffset + 0 * cont2csize + 26] = current_source[soffset + 728];
          current_target[toffset + 1 * cont2csize + 26] = current_source[soffset + 729];
          current_target[toffset + 2 * cont2csize + 26] = current_source[soffset + 730];
          current_target[toffset + 3 * cont2csize + 26] = current_source[soffset + 731];
          current_target[toffset + 4 * cont2csize + 26] = current_source[soffset + 732];
          current_target[toffset + 5 * cont2csize + 26] = current_source[soffset + 733];
          current_target[toffset + 6 * cont2csize + 26] = current_source[soffset + 734];
          current_target[toffset + 7 * cont2csize + 26] = current_source[soffset + 735];
          current_target[toffset + 8 * cont2csize + 26] = current_source[soffset + 736];
          current_target[toffset + 9 * cont2csize + 26] = current_source[soffset + 737];
          current_target[toffset + 10 * cont2csize + 26] = current_source[soffset + 738];
          current_target[toffset + 11 * cont2csize + 26] = current_source[soffset + 739];
          current_target[toffset + 12 * cont2csize + 26] = current_source[soffset + 740];
          current_target[toffset + 13 * cont2csize + 26] = current_source[soffset + 741];
          current_target[toffset + 14 * cont2csize + 26] = current_source[soffset + 742];
          current_target[toffset + 15 * cont2csize + 26] = current_source[soffset + 743];
          current_target[toffset + 16 * cont2csize + 26] = current_source[soffset + 744];
          current_target[toffset + 17 * cont2csize + 26] = current_source[soffset + 745];
          current_target[toffset + 18 * cont2csize + 26] = current_source[soffset + 746];
          current_target[toffset + 19 * cont2csize + 26] = current_source[soffset + 747];
          current_target[toffset + 20 * cont2csize + 26] = current_source[soffset + 748];
          current_target[toffset + 21 * cont2csize + 26] = current_source[soffset + 749];
          current_target[toffset + 22 * cont2csize + 26] = current_source[soffset + 750];
          current_target[toffset + 23 * cont2csize + 26] = current_source[soffset + 751];
          current_target[toffset + 24 * cont2csize + 26] = current_source[soffset + 752];
          current_target[toffset + 25 * cont2csize + 26] = current_source[soffset + 753];
          current_target[toffset + 26 * cont2csize + 26] = current_source[soffset + 754];
          current_target[toffset + 27 * cont2csize + 26] = current_source[soffset + 755];
          current_target[toffset + 0 * cont2csize + 27] = current_source[soffset + 756];
          current_target[toffset + 1 * cont2csize + 27] = current_source[soffset + 757];
          current_target[toffset + 2 * cont2csize + 27] = current_source[soffset + 758];
          current_target[toffset + 3 * cont2csize + 27] = current_source[soffset + 759];
          current_target[toffset + 4 * cont2csize + 27] = current_source[soffset + 760];
          current_target[toffset + 5 * cont2csize + 27] = current_source[soffset + 761];
          current_target[toffset + 6 * cont2csize + 27] = current_source[soffset + 762];
          current_target[toffset + 7 * cont2csize + 27] = current_source[soffset + 763];
          current_target[toffset + 8 * cont2csize + 27] = current_source[soffset + 764];
          current_target[toffset + 9 * cont2csize + 27] = current_source[soffset + 765];
          current_target[toffset + 10 * cont2csize + 27] = current_source[soffset + 766];
          current_target[toffset + 11 * cont2csize + 27] = current_source[soffset + 767];
          current_target[toffset + 12 * cont2csize + 27] = current_source[soffset + 768];
          current_target[toffset + 13 * cont2csize + 27] = current_source[soffset + 769];
          current_target[toffset + 14 * cont2csize + 27] = current_source[soffset + 770];
          current_target[toffset + 15 * cont2csize + 27] = current_source[soffset + 771];
          current_target[toffset + 16 * cont2csize + 27] = current_source[soffset + 772];
          current_target[toffset + 17 * cont2csize + 27] = current_source[soffset + 773];
          current_target[toffset + 18 * cont2csize + 27] = current_source[soffset + 774];
          current_target[toffset + 19 * cont2csize + 27] = current_source[soffset + 775];
          current_target[toffset + 20 * cont2csize + 27] = current_source[soffset + 776];
          current_target[toffset + 21 * cont2csize + 27] = current_source[soffset + 777];
          current_target[toffset + 22 * cont2csize + 27] = current_source[soffset + 778];
          current_target[toffset + 23 * cont2csize + 27] = current_source[soffset + 779];
          current_target[toffset + 24 * cont2csize + 27] = current_source[soffset + 780];
          current_target[toffset + 25 * cont2csize + 27] = current_source[soffset + 781];
          current_target[toffset + 26 * cont2csize + 27] = current_source[soffset + 782];
          current_target[toffset + 27 * cont2csize + 27] = current_source[soffset + 783];
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = 28 * c3end;
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      complex<double>* current_target = &target[offset];
      const complex<double>* current_source = &source[offset];

      for (int c2 = 0; c2 != c2end; ++c2) {
        for (int c3 = 0; c3 != c3end; ++c3) {
          const int c3x3end = c3 * 28;
          const int soffset = 784 * (c3 + c3end * c2);
          const int toffset = 28 * c2 * cont3csize + c3x3end;
          copy_n(current_source+soffset+  0,  28, current_target+toffset+ 0*cont3csize);
          copy_n(current_source+soffset+ 28,  28, current_target+toffset+ 1*cont3csize);
          copy_n(current_source+soffset+ 56,  28, current_target+toffset+ 2*cont3csize);
          copy_n(current_source+soffset+ 84,  28, current_target+toffset+ 3*cont3csize);
          copy_n(current_source+soffset+112,  28, current_target+toffset+ 4*cont3csize);
          copy_n(current_source+soffset+140,  28, current_target+toffset+ 5*cont3csize);
          copy_n(current_source+soffset+168,  28, current_target+toffset+ 6*cont3csize);
          copy_n(current_source+soffset+196,  28, current_target+toffset+ 7*cont3csize);
          copy_n(current_source+soffset+224,  28, current_target+toffset+ 8*cont3csize);
          copy_n(current_source+soffset+252,  28, current_target+toffset+ 9*cont3csize);
          copy_n(current_source+soffset+280,  28, current_target+toffset+10*cont3csize);
          copy_n(current_source+soffset+308,  28, current_target+toffset+11*cont3csize);
          copy_n(current_source+soffset+336,  28, current_target+toffset+12*cont3csize);
          copy_n(current_source+soffset+364,  28, current_target+toffset+13*cont3csize);
          copy_n(current_source+soffset+392,  28, current_target+toffset+14*cont3csize);
          copy_n(current_source+soffset+420,  28, current_target+toffset+15*cont3csize);
          copy_n(current_source+soffset+448,  28, current_target+toffset+16*cont3csize);
          copy_n(current_source+soffset+476,  28, current_target+toffset+17*cont3csize);
          copy_n(current_source+soffset+504,  28, current_target+toffset+18*cont3csize);
          copy_n(current_source+soffset+532,  28, current_target+toffset+19*cont3csize);
          copy_n(current_source+soffset+560,  28, current_target+toffset+20*cont3csize);
          copy_n(current_source+soffset+588,  28, current_target+toffset+21*cont3csize);
          copy_n(current_source+soffset+616,  28, current_target+toffset+22*cont3csize);
          copy_n(current_source+soffset+644,  28, current_target+toffset+23*cont3csize);
          copy_n(current_source+soffset+672,  28, current_target+toffset+24*cont3csize);
          copy_n(current_source+soffset+700,  28, current_target+toffset+25*cont3csize);
          copy_n(current_source+soffset+728,  28, current_target+toffset+26*cont3csize);
          copy_n(current_source+soffset+756,  28, current_target+toffset+27*cont3csize);
        }
      }

    }
  }

}



