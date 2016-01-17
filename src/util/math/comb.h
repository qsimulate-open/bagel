//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: comb.h
// Copyright (C) 2011 Toru Shiozaki
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

//

#ifndef __BAGEL_UTIL_COMB_H
#define __BAGEL_UTIL_COMB_H

#include <algorithm>
#include <array>
#include <cassert>

namespace bagel {

class Comb {
    constexpr static int max_ = 65;
    std::array<size_t, max_*max_> c_;
  public:
    Comb() {
      std::fill(c_.begin(), c_.end(), 0);
      c_[ 0*max_+ 0] =                    1ull; c_[ 1*max_+ 0] =                    1ull; c_[ 1*max_+ 1] =                    1ull; c_[ 2*max_+ 0] =                    1ull;
      c_[ 2*max_+ 1] =                    2ull; c_[ 2*max_+ 2] =                    1ull; c_[ 3*max_+ 0] =                    1ull; c_[ 3*max_+ 1] =                    3ull;
      c_[ 3*max_+ 2] =                    3ull; c_[ 3*max_+ 3] =                    1ull; c_[ 4*max_+ 0] =                    1ull; c_[ 4*max_+ 1] =                    4ull;
      c_[ 4*max_+ 2] =                    6ull; c_[ 4*max_+ 3] =                    4ull; c_[ 4*max_+ 4] =                    1ull; c_[ 5*max_+ 0] =                    1ull;
      c_[ 5*max_+ 1] =                    5ull; c_[ 5*max_+ 2] =                   10ull; c_[ 5*max_+ 3] =                   10ull; c_[ 5*max_+ 4] =                    5ull;
      c_[ 5*max_+ 5] =                    1ull; c_[ 6*max_+ 0] =                    1ull; c_[ 6*max_+ 1] =                    6ull; c_[ 6*max_+ 2] =                   15ull;
      c_[ 6*max_+ 3] =                   20ull; c_[ 6*max_+ 4] =                   15ull; c_[ 6*max_+ 5] =                    6ull; c_[ 6*max_+ 6] =                    1ull;
      c_[ 7*max_+ 0] =                    1ull; c_[ 7*max_+ 1] =                    7ull; c_[ 7*max_+ 2] =                   21ull; c_[ 7*max_+ 3] =                   35ull;
      c_[ 7*max_+ 4] =                   35ull; c_[ 7*max_+ 5] =                   21ull; c_[ 7*max_+ 6] =                    7ull; c_[ 7*max_+ 7] =                    1ull;
      c_[ 8*max_+ 0] =                    1ull; c_[ 8*max_+ 1] =                    8ull; c_[ 8*max_+ 2] =                   28ull; c_[ 8*max_+ 3] =                   56ull;
      c_[ 8*max_+ 4] =                   70ull; c_[ 8*max_+ 5] =                   56ull; c_[ 8*max_+ 6] =                   28ull; c_[ 8*max_+ 7] =                    8ull;
      c_[ 8*max_+ 8] =                    1ull; c_[ 9*max_+ 0] =                    1ull; c_[ 9*max_+ 1] =                    9ull; c_[ 9*max_+ 2] =                   36ull;
      c_[ 9*max_+ 3] =                   84ull; c_[ 9*max_+ 4] =                  126ull; c_[ 9*max_+ 5] =                  126ull; c_[ 9*max_+ 6] =                   84ull;
      c_[ 9*max_+ 7] =                   36ull; c_[ 9*max_+ 8] =                    9ull; c_[ 9*max_+ 9] =                    1ull; c_[10*max_+ 0] =                    1ull;
      c_[10*max_+ 1] =                   10ull; c_[10*max_+ 2] =                   45ull; c_[10*max_+ 3] =                  120ull; c_[10*max_+ 4] =                  210ull;
      c_[10*max_+ 5] =                  252ull; c_[10*max_+ 6] =                  210ull; c_[10*max_+ 7] =                  120ull; c_[10*max_+ 8] =                   45ull;
      c_[10*max_+ 9] =                   10ull; c_[10*max_+10] =                    1ull; c_[11*max_+ 0] =                    1ull; c_[11*max_+ 1] =                   11ull;
      c_[11*max_+ 2] =                   55ull; c_[11*max_+ 3] =                  165ull; c_[11*max_+ 4] =                  330ull; c_[11*max_+ 5] =                  462ull;
      c_[11*max_+ 6] =                  462ull; c_[11*max_+ 7] =                  330ull; c_[11*max_+ 8] =                  165ull; c_[11*max_+ 9] =                   55ull;
      c_[11*max_+10] =                   11ull; c_[11*max_+11] =                    1ull; c_[12*max_+ 0] =                    1ull; c_[12*max_+ 1] =                   12ull;
      c_[12*max_+ 2] =                   66ull; c_[12*max_+ 3] =                  220ull; c_[12*max_+ 4] =                  495ull; c_[12*max_+ 5] =                  792ull;
      c_[12*max_+ 6] =                  924ull; c_[12*max_+ 7] =                  792ull; c_[12*max_+ 8] =                  495ull; c_[12*max_+ 9] =                  220ull;
      c_[12*max_+10] =                   66ull; c_[12*max_+11] =                   12ull; c_[12*max_+12] =                    1ull; c_[13*max_+ 0] =                    1ull;
      c_[13*max_+ 1] =                   13ull; c_[13*max_+ 2] =                   78ull; c_[13*max_+ 3] =                  286ull; c_[13*max_+ 4] =                  715ull;
      c_[13*max_+ 5] =                 1287ull; c_[13*max_+ 6] =                 1716ull; c_[13*max_+ 7] =                 1716ull; c_[13*max_+ 8] =                 1287ull;
      c_[13*max_+ 9] =                  715ull; c_[13*max_+10] =                  286ull; c_[13*max_+11] =                   78ull; c_[13*max_+12] =                   13ull;
      c_[13*max_+13] =                    1ull; c_[14*max_+ 0] =                    1ull; c_[14*max_+ 1] =                   14ull; c_[14*max_+ 2] =                   91ull;
      c_[14*max_+ 3] =                  364ull; c_[14*max_+ 4] =                 1001ull; c_[14*max_+ 5] =                 2002ull; c_[14*max_+ 6] =                 3003ull;
      c_[14*max_+ 7] =                 3432ull; c_[14*max_+ 8] =                 3003ull; c_[14*max_+ 9] =                 2002ull; c_[14*max_+10] =                 1001ull;
      c_[14*max_+11] =                  364ull; c_[14*max_+12] =                   91ull; c_[14*max_+13] =                   14ull; c_[14*max_+14] =                    1ull;
      c_[15*max_+ 0] =                    1ull; c_[15*max_+ 1] =                   15ull; c_[15*max_+ 2] =                  105ull; c_[15*max_+ 3] =                  455ull;
      c_[15*max_+ 4] =                 1365ull; c_[15*max_+ 5] =                 3003ull; c_[15*max_+ 6] =                 5005ull; c_[15*max_+ 7] =                 6435ull;
      c_[15*max_+ 8] =                 6435ull; c_[15*max_+ 9] =                 5005ull; c_[15*max_+10] =                 3003ull; c_[15*max_+11] =                 1365ull;
      c_[15*max_+12] =                  455ull; c_[15*max_+13] =                  105ull; c_[15*max_+14] =                   15ull; c_[15*max_+15] =                    1ull;
      c_[16*max_+ 0] =                    1ull; c_[16*max_+ 1] =                   16ull; c_[16*max_+ 2] =                  120ull; c_[16*max_+ 3] =                  560ull;
      c_[16*max_+ 4] =                 1820ull; c_[16*max_+ 5] =                 4368ull; c_[16*max_+ 6] =                 8008ull; c_[16*max_+ 7] =                11440ull;
      c_[16*max_+ 8] =                12870ull; c_[16*max_+ 9] =                11440ull; c_[16*max_+10] =                 8008ull; c_[16*max_+11] =                 4368ull;
      c_[16*max_+12] =                 1820ull; c_[16*max_+13] =                  560ull; c_[16*max_+14] =                  120ull; c_[16*max_+15] =                   16ull;
      c_[16*max_+16] =                    1ull; c_[17*max_+ 0] =                    1ull; c_[17*max_+ 1] =                   17ull; c_[17*max_+ 2] =                  136ull;
      c_[17*max_+ 3] =                  680ull; c_[17*max_+ 4] =                 2380ull; c_[17*max_+ 5] =                 6188ull; c_[17*max_+ 6] =                12376ull;
      c_[17*max_+ 7] =                19448ull; c_[17*max_+ 8] =                24310ull; c_[17*max_+ 9] =                24310ull; c_[17*max_+10] =                19448ull;
      c_[17*max_+11] =                12376ull; c_[17*max_+12] =                 6188ull; c_[17*max_+13] =                 2380ull; c_[17*max_+14] =                  680ull;
      c_[17*max_+15] =                  136ull; c_[17*max_+16] =                   17ull; c_[17*max_+17] =                    1ull; c_[18*max_+ 0] =                    1ull;
      c_[18*max_+ 1] =                   18ull; c_[18*max_+ 2] =                  153ull; c_[18*max_+ 3] =                  816ull; c_[18*max_+ 4] =                 3060ull;
      c_[18*max_+ 5] =                 8568ull; c_[18*max_+ 6] =                18564ull; c_[18*max_+ 7] =                31824ull; c_[18*max_+ 8] =                43758ull;
      c_[18*max_+ 9] =                48620ull; c_[18*max_+10] =                43758ull; c_[18*max_+11] =                31824ull; c_[18*max_+12] =                18564ull;
      c_[18*max_+13] =                 8568ull; c_[18*max_+14] =                 3060ull; c_[18*max_+15] =                  816ull; c_[18*max_+16] =                  153ull;
      c_[18*max_+17] =                   18ull; c_[18*max_+18] =                    1ull; c_[19*max_+ 0] =                    1ull; c_[19*max_+ 1] =                   19ull;
      c_[19*max_+ 2] =                  171ull; c_[19*max_+ 3] =                  969ull; c_[19*max_+ 4] =                 3876ull; c_[19*max_+ 5] =                11628ull;
      c_[19*max_+ 6] =                27132ull; c_[19*max_+ 7] =                50388ull; c_[19*max_+ 8] =                75582ull; c_[19*max_+ 9] =                92378ull;
      c_[19*max_+10] =                92378ull; c_[19*max_+11] =                75582ull; c_[19*max_+12] =                50388ull; c_[19*max_+13] =                27132ull;
      c_[19*max_+14] =                11628ull; c_[19*max_+15] =                 3876ull; c_[19*max_+16] =                  969ull; c_[19*max_+17] =                  171ull;
      c_[19*max_+18] =                   19ull; c_[19*max_+19] =                    1ull; c_[20*max_+ 0] =                    1ull; c_[20*max_+ 1] =                   20ull;
      c_[20*max_+ 2] =                  190ull; c_[20*max_+ 3] =                 1140ull; c_[20*max_+ 4] =                 4845ull; c_[20*max_+ 5] =                15504ull;
      c_[20*max_+ 6] =                38760ull; c_[20*max_+ 7] =                77520ull; c_[20*max_+ 8] =               125970ull; c_[20*max_+ 9] =               167960ull;
      c_[20*max_+10] =               184756ull; c_[20*max_+11] =               167960ull; c_[20*max_+12] =               125970ull; c_[20*max_+13] =                77520ull;
      c_[20*max_+14] =                38760ull; c_[20*max_+15] =                15504ull; c_[20*max_+16] =                 4845ull; c_[20*max_+17] =                 1140ull;
      c_[20*max_+18] =                  190ull; c_[20*max_+19] =                   20ull; c_[20*max_+20] =                    1ull; c_[21*max_+ 0] =                    1ull;
      c_[21*max_+ 1] =                   21ull; c_[21*max_+ 2] =                  210ull; c_[21*max_+ 3] =                 1330ull; c_[21*max_+ 4] =                 5985ull;
      c_[21*max_+ 5] =                20349ull; c_[21*max_+ 6] =                54264ull; c_[21*max_+ 7] =               116280ull; c_[21*max_+ 8] =               203490ull;
      c_[21*max_+ 9] =               293930ull; c_[21*max_+10] =               352716ull; c_[21*max_+11] =               352716ull; c_[21*max_+12] =               293930ull;
      c_[21*max_+13] =               203490ull; c_[21*max_+14] =               116280ull; c_[21*max_+15] =                54264ull; c_[21*max_+16] =                20349ull;
      c_[21*max_+17] =                 5985ull; c_[21*max_+18] =                 1330ull; c_[21*max_+19] =                  210ull; c_[21*max_+20] =                   21ull;
      c_[21*max_+21] =                    1ull; c_[22*max_+ 0] =                    1ull; c_[22*max_+ 1] =                   22ull; c_[22*max_+ 2] =                  231ull;
      c_[22*max_+ 3] =                 1540ull; c_[22*max_+ 4] =                 7315ull; c_[22*max_+ 5] =                26334ull; c_[22*max_+ 6] =                74613ull;
      c_[22*max_+ 7] =               170544ull; c_[22*max_+ 8] =               319770ull; c_[22*max_+ 9] =               497420ull; c_[22*max_+10] =               646646ull;
      c_[22*max_+11] =               705432ull; c_[22*max_+12] =               646646ull; c_[22*max_+13] =               497420ull; c_[22*max_+14] =               319770ull;
      c_[22*max_+15] =               170544ull; c_[22*max_+16] =                74613ull; c_[22*max_+17] =                26334ull; c_[22*max_+18] =                 7315ull;
      c_[22*max_+19] =                 1540ull; c_[22*max_+20] =                  231ull; c_[22*max_+21] =                   22ull; c_[22*max_+22] =                    1ull;
      c_[23*max_+ 0] =                    1ull; c_[23*max_+ 1] =                   23ull; c_[23*max_+ 2] =                  253ull; c_[23*max_+ 3] =                 1771ull;
      c_[23*max_+ 4] =                 8855ull; c_[23*max_+ 5] =                33649ull; c_[23*max_+ 6] =               100947ull; c_[23*max_+ 7] =               245157ull;
      c_[23*max_+ 8] =               490314ull; c_[23*max_+ 9] =               817190ull; c_[23*max_+10] =              1144066ull; c_[23*max_+11] =              1352078ull;
      c_[23*max_+12] =              1352078ull; c_[23*max_+13] =              1144066ull; c_[23*max_+14] =               817190ull; c_[23*max_+15] =               490314ull;
      c_[23*max_+16] =               245157ull; c_[23*max_+17] =               100947ull; c_[23*max_+18] =                33649ull; c_[23*max_+19] =                 8855ull;
      c_[23*max_+20] =                 1771ull; c_[23*max_+21] =                  253ull; c_[23*max_+22] =                   23ull; c_[23*max_+23] =                    1ull;
      c_[24*max_+ 0] =                    1ull; c_[24*max_+ 1] =                   24ull; c_[24*max_+ 2] =                  276ull; c_[24*max_+ 3] =                 2024ull;
      c_[24*max_+ 4] =                10626ull; c_[24*max_+ 5] =                42504ull; c_[24*max_+ 6] =               134596ull; c_[24*max_+ 7] =               346104ull;
      c_[24*max_+ 8] =               735471ull; c_[24*max_+ 9] =              1307504ull; c_[24*max_+10] =              1961256ull; c_[24*max_+11] =              2496144ull;
      c_[24*max_+12] =              2704156ull; c_[24*max_+13] =              2496144ull; c_[24*max_+14] =              1961256ull; c_[24*max_+15] =              1307504ull;
      c_[24*max_+16] =               735471ull; c_[24*max_+17] =               346104ull; c_[24*max_+18] =               134596ull; c_[24*max_+19] =                42504ull;
      c_[24*max_+20] =                10626ull; c_[24*max_+21] =                 2024ull; c_[24*max_+22] =                  276ull; c_[24*max_+23] =                   24ull;
      c_[24*max_+24] =                    1ull; c_[25*max_+ 0] =                    1ull; c_[25*max_+ 1] =                   25ull; c_[25*max_+ 2] =                  300ull;
      c_[25*max_+ 3] =                 2300ull; c_[25*max_+ 4] =                12650ull; c_[25*max_+ 5] =                53130ull; c_[25*max_+ 6] =               177100ull;
      c_[25*max_+ 7] =               480700ull; c_[25*max_+ 8] =              1081575ull; c_[25*max_+ 9] =              2042975ull; c_[25*max_+10] =              3268760ull;
      c_[25*max_+11] =              4457400ull; c_[25*max_+12] =              5200300ull; c_[25*max_+13] =              5200300ull; c_[25*max_+14] =              4457400ull;
      c_[25*max_+15] =              3268760ull; c_[25*max_+16] =              2042975ull; c_[25*max_+17] =              1081575ull; c_[25*max_+18] =               480700ull;
      c_[25*max_+19] =               177100ull; c_[25*max_+20] =                53130ull; c_[25*max_+21] =                12650ull; c_[25*max_+22] =                 2300ull;
      c_[25*max_+23] =                  300ull; c_[25*max_+24] =                   25ull; c_[25*max_+25] =                    1ull; c_[26*max_+ 0] =                    1ull;
      c_[26*max_+ 1] =                   26ull; c_[26*max_+ 2] =                  325ull; c_[26*max_+ 3] =                 2600ull; c_[26*max_+ 4] =                14950ull;
      c_[26*max_+ 5] =                65780ull; c_[26*max_+ 6] =               230230ull; c_[26*max_+ 7] =               657800ull; c_[26*max_+ 8] =              1562275ull;
      c_[26*max_+ 9] =              3124550ull; c_[26*max_+10] =              5311735ull; c_[26*max_+11] =              7726160ull; c_[26*max_+12] =              9657700ull;
      c_[26*max_+13] =             10400600ull; c_[26*max_+14] =              9657700ull; c_[26*max_+15] =              7726160ull; c_[26*max_+16] =              5311735ull;
      c_[26*max_+17] =              3124550ull; c_[26*max_+18] =              1562275ull; c_[26*max_+19] =               657800ull; c_[26*max_+20] =               230230ull;
      c_[26*max_+21] =                65780ull; c_[26*max_+22] =                14950ull; c_[26*max_+23] =                 2600ull; c_[26*max_+24] =                  325ull;
      c_[26*max_+25] =                   26ull; c_[26*max_+26] =                    1ull; c_[27*max_+ 0] =                    1ull; c_[27*max_+ 1] =                   27ull;
      c_[27*max_+ 2] =                  351ull; c_[27*max_+ 3] =                 2925ull; c_[27*max_+ 4] =                17550ull; c_[27*max_+ 5] =                80730ull;
      c_[27*max_+ 6] =               296010ull; c_[27*max_+ 7] =               888030ull; c_[27*max_+ 8] =              2220075ull; c_[27*max_+ 9] =              4686825ull;
      c_[27*max_+10] =              8436285ull; c_[27*max_+11] =             13037895ull; c_[27*max_+12] =             17383860ull; c_[27*max_+13] =             20058300ull;
      c_[27*max_+14] =             20058300ull; c_[27*max_+15] =             17383860ull; c_[27*max_+16] =             13037895ull; c_[27*max_+17] =              8436285ull;
      c_[27*max_+18] =              4686825ull; c_[27*max_+19] =              2220075ull; c_[27*max_+20] =               888030ull; c_[27*max_+21] =               296010ull;
      c_[27*max_+22] =                80730ull; c_[27*max_+23] =                17550ull; c_[27*max_+24] =                 2925ull; c_[27*max_+25] =                  351ull;
      c_[27*max_+26] =                   27ull; c_[27*max_+27] =                    1ull; c_[28*max_+ 0] =                    1ull; c_[28*max_+ 1] =                   28ull;
      c_[28*max_+ 2] =                  378ull; c_[28*max_+ 3] =                 3276ull; c_[28*max_+ 4] =                20475ull; c_[28*max_+ 5] =                98280ull;
      c_[28*max_+ 6] =               376740ull; c_[28*max_+ 7] =              1184040ull; c_[28*max_+ 8] =              3108105ull; c_[28*max_+ 9] =              6906900ull;
      c_[28*max_+10] =             13123110ull; c_[28*max_+11] =             21474180ull; c_[28*max_+12] =             30421755ull; c_[28*max_+13] =             37442160ull;
      c_[28*max_+14] =             40116600ull; c_[28*max_+15] =             37442160ull; c_[28*max_+16] =             30421755ull; c_[28*max_+17] =             21474180ull;
      c_[28*max_+18] =             13123110ull; c_[28*max_+19] =              6906900ull; c_[28*max_+20] =              3108105ull; c_[28*max_+21] =              1184040ull;
      c_[28*max_+22] =               376740ull; c_[28*max_+23] =                98280ull; c_[28*max_+24] =                20475ull; c_[28*max_+25] =                 3276ull;
      c_[28*max_+26] =                  378ull; c_[28*max_+27] =                   28ull; c_[28*max_+28] =                    1ull; c_[29*max_+ 0] =                    1ull;
      c_[29*max_+ 1] =                   29ull; c_[29*max_+ 2] =                  406ull; c_[29*max_+ 3] =                 3654ull; c_[29*max_+ 4] =                23751ull;
      c_[29*max_+ 5] =               118755ull; c_[29*max_+ 6] =               475020ull; c_[29*max_+ 7] =              1560780ull; c_[29*max_+ 8] =              4292145ull;
      c_[29*max_+ 9] =             10015005ull; c_[29*max_+10] =             20030010ull; c_[29*max_+11] =             34597290ull; c_[29*max_+12] =             51895935ull;
      c_[29*max_+13] =             67863915ull; c_[29*max_+14] =             77558760ull; c_[29*max_+15] =             77558760ull; c_[29*max_+16] =             67863915ull;
      c_[29*max_+17] =             51895935ull; c_[29*max_+18] =             34597290ull; c_[29*max_+19] =             20030010ull; c_[29*max_+20] =             10015005ull;
      c_[29*max_+21] =              4292145ull; c_[29*max_+22] =              1560780ull; c_[29*max_+23] =               475020ull; c_[29*max_+24] =               118755ull;
      c_[29*max_+25] =                23751ull; c_[29*max_+26] =                 3654ull; c_[29*max_+27] =                  406ull; c_[29*max_+28] =                   29ull;
      c_[29*max_+29] =                    1ull; c_[30*max_+ 0] =                    1ull; c_[30*max_+ 1] =                   30ull; c_[30*max_+ 2] =                  435ull;
      c_[30*max_+ 3] =                 4060ull; c_[30*max_+ 4] =                27405ull; c_[30*max_+ 5] =               142506ull; c_[30*max_+ 6] =               593775ull;
      c_[30*max_+ 7] =              2035800ull; c_[30*max_+ 8] =              5852925ull; c_[30*max_+ 9] =             14307150ull; c_[30*max_+10] =             30045015ull;
      c_[30*max_+11] =             54627300ull; c_[30*max_+12] =             86493225ull; c_[30*max_+13] =            119759850ull; c_[30*max_+14] =            145422675ull;
      c_[30*max_+15] =            155117520ull; c_[30*max_+16] =            145422675ull; c_[30*max_+17] =            119759850ull; c_[30*max_+18] =             86493225ull;
      c_[30*max_+19] =             54627300ull; c_[30*max_+20] =             30045015ull; c_[30*max_+21] =             14307150ull; c_[30*max_+22] =              5852925ull;
      c_[30*max_+23] =              2035800ull; c_[30*max_+24] =               593775ull; c_[30*max_+25] =               142506ull; c_[30*max_+26] =                27405ull;
      c_[30*max_+27] =                 4060ull; c_[30*max_+28] =                  435ull; c_[30*max_+29] =                   30ull; c_[30*max_+30] =                    1ull;
      c_[31*max_+ 0] =                    1ull; c_[31*max_+ 1] =                   31ull; c_[31*max_+ 2] =                  465ull; c_[31*max_+ 3] =                 4495ull;
      c_[31*max_+ 4] =                31465ull; c_[31*max_+ 5] =               169911ull; c_[31*max_+ 6] =               736281ull; c_[31*max_+ 7] =              2629575ull;
      c_[31*max_+ 8] =              7888725ull; c_[31*max_+ 9] =             20160075ull; c_[31*max_+10] =             44352165ull; c_[31*max_+11] =             84672315ull;
      c_[31*max_+12] =            141120525ull; c_[31*max_+13] =            206253075ull; c_[31*max_+14] =            265182525ull; c_[31*max_+15] =            300540195ull;
      c_[31*max_+16] =            300540195ull; c_[31*max_+17] =            265182525ull; c_[31*max_+18] =            206253075ull; c_[31*max_+19] =            141120525ull;
      c_[31*max_+20] =             84672315ull; c_[31*max_+21] =             44352165ull; c_[31*max_+22] =             20160075ull; c_[31*max_+23] =              7888725ull;
      c_[31*max_+24] =              2629575ull; c_[31*max_+25] =               736281ull; c_[31*max_+26] =               169911ull; c_[31*max_+27] =                31465ull;
      c_[31*max_+28] =                 4495ull; c_[31*max_+29] =                  465ull; c_[31*max_+30] =                   31ull; c_[31*max_+31] =                    1ull;
      c_[32*max_+ 0] =                    1ull; c_[32*max_+ 1] =                   32ull; c_[32*max_+ 2] =                  496ull; c_[32*max_+ 3] =                 4960ull;
      c_[32*max_+ 4] =                35960ull; c_[32*max_+ 5] =               201376ull; c_[32*max_+ 6] =               906192ull; c_[32*max_+ 7] =              3365856ull;
      c_[32*max_+ 8] =             10518300ull; c_[32*max_+ 9] =             28048800ull; c_[32*max_+10] =             64512240ull; c_[32*max_+11] =            129024480ull;
      c_[32*max_+12] =            225792840ull; c_[32*max_+13] =            347373600ull; c_[32*max_+14] =            471435600ull; c_[32*max_+15] =            565722720ull;
      c_[32*max_+16] =            601080390ull; c_[32*max_+17] =            565722720ull; c_[32*max_+18] =            471435600ull; c_[32*max_+19] =            347373600ull;
      c_[32*max_+20] =            225792840ull; c_[32*max_+21] =            129024480ull; c_[32*max_+22] =             64512240ull; c_[32*max_+23] =             28048800ull;
      c_[32*max_+24] =             10518300ull; c_[32*max_+25] =              3365856ull; c_[32*max_+26] =               906192ull; c_[32*max_+27] =               201376ull;
      c_[32*max_+28] =                35960ull; c_[32*max_+29] =                 4960ull; c_[32*max_+30] =                  496ull; c_[32*max_+31] =                   32ull;
      c_[32*max_+32] =                    1ull; c_[33*max_+ 0] =                    1ull; c_[33*max_+ 1] =                   33ull; c_[33*max_+ 2] =                  528ull;
      c_[33*max_+ 3] =                 5456ull; c_[33*max_+ 4] =                40920ull; c_[33*max_+ 5] =               237336ull; c_[33*max_+ 6] =              1107568ull;
      c_[33*max_+ 7] =              4272048ull; c_[33*max_+ 8] =             13884156ull; c_[33*max_+ 9] =             38567100ull; c_[33*max_+10] =             92561040ull;
      c_[33*max_+11] =            193536720ull; c_[33*max_+12] =            354817320ull; c_[33*max_+13] =            573166440ull; c_[33*max_+14] =            818809200ull;
      c_[33*max_+15] =           1037158320ull; c_[33*max_+16] =           1166803110ull; c_[33*max_+17] =           1166803110ull; c_[33*max_+18] =           1037158320ull;
      c_[33*max_+19] =            818809200ull; c_[33*max_+20] =            573166440ull; c_[33*max_+21] =            354817320ull; c_[33*max_+22] =            193536720ull;
      c_[33*max_+23] =             92561040ull; c_[33*max_+24] =             38567100ull; c_[33*max_+25] =             13884156ull; c_[33*max_+26] =              4272048ull;
      c_[33*max_+27] =              1107568ull; c_[33*max_+28] =               237336ull; c_[33*max_+29] =                40920ull; c_[33*max_+30] =                 5456ull;
      c_[33*max_+31] =                  528ull; c_[33*max_+32] =                   33ull; c_[33*max_+33] =                    1ull; c_[34*max_+ 0] =                    1ull;
      c_[34*max_+ 1] =                   34ull; c_[34*max_+ 2] =                  561ull; c_[34*max_+ 3] =                 5984ull; c_[34*max_+ 4] =                46376ull;
      c_[34*max_+ 5] =               278256ull; c_[34*max_+ 6] =              1344904ull; c_[34*max_+ 7] =              5379616ull; c_[34*max_+ 8] =             18156204ull;
      c_[34*max_+ 9] =             52451256ull; c_[34*max_+10] =            131128140ull; c_[34*max_+11] =            286097760ull; c_[34*max_+12] =            548354040ull;
      c_[34*max_+13] =            927983760ull; c_[34*max_+14] =           1391975640ull; c_[34*max_+15] =           1855967520ull; c_[34*max_+16] =           2203961430ull;
      c_[34*max_+17] =           2333606220ull; c_[34*max_+18] =           2203961430ull; c_[34*max_+19] =           1855967520ull; c_[34*max_+20] =           1391975640ull;
      c_[34*max_+21] =            927983760ull; c_[34*max_+22] =            548354040ull; c_[34*max_+23] =            286097760ull; c_[34*max_+24] =            131128140ull;
      c_[34*max_+25] =             52451256ull; c_[34*max_+26] =             18156204ull; c_[34*max_+27] =              5379616ull; c_[34*max_+28] =              1344904ull;
      c_[34*max_+29] =               278256ull; c_[34*max_+30] =                46376ull; c_[34*max_+31] =                 5984ull; c_[34*max_+32] =                  561ull;
      c_[34*max_+33] =                   34ull; c_[34*max_+34] =                    1ull; c_[35*max_+ 0] =                    1ull; c_[35*max_+ 1] =                   35ull;
      c_[35*max_+ 2] =                  595ull; c_[35*max_+ 3] =                 6545ull; c_[35*max_+ 4] =                52360ull; c_[35*max_+ 5] =               324632ull;
      c_[35*max_+ 6] =              1623160ull; c_[35*max_+ 7] =              6724520ull; c_[35*max_+ 8] =             23535820ull; c_[35*max_+ 9] =             70607460ull;
      c_[35*max_+10] =            183579396ull; c_[35*max_+11] =            417225900ull; c_[35*max_+12] =            834451800ull; c_[35*max_+13] =           1476337800ull;
      c_[35*max_+14] =           2319959400ull; c_[35*max_+15] =           3247943160ull; c_[35*max_+16] =           4059928950ull; c_[35*max_+17] =           4537567650ull;
      c_[35*max_+18] =           4537567650ull; c_[35*max_+19] =           4059928950ull; c_[35*max_+20] =           3247943160ull; c_[35*max_+21] =           2319959400ull;
      c_[35*max_+22] =           1476337800ull; c_[35*max_+23] =            834451800ull; c_[35*max_+24] =            417225900ull; c_[35*max_+25] =            183579396ull;
      c_[35*max_+26] =             70607460ull; c_[35*max_+27] =             23535820ull; c_[35*max_+28] =              6724520ull; c_[35*max_+29] =              1623160ull;
      c_[35*max_+30] =               324632ull; c_[35*max_+31] =                52360ull; c_[35*max_+32] =                 6545ull; c_[35*max_+33] =                  595ull;
      c_[35*max_+34] =                   35ull; c_[35*max_+35] =                    1ull; c_[36*max_+ 0] =                    1ull; c_[36*max_+ 1] =                   36ull;
      c_[36*max_+ 2] =                  630ull; c_[36*max_+ 3] =                 7140ull; c_[36*max_+ 4] =                58905ull; c_[36*max_+ 5] =               376992ull;
      c_[36*max_+ 6] =              1947792ull; c_[36*max_+ 7] =              8347680ull; c_[36*max_+ 8] =             30260340ull; c_[36*max_+ 9] =             94143280ull;
      c_[36*max_+10] =            254186856ull; c_[36*max_+11] =            600805296ull; c_[36*max_+12] =           1251677700ull; c_[36*max_+13] =           2310789600ull;
      c_[36*max_+14] =           3796297200ull; c_[36*max_+15] =           5567902560ull; c_[36*max_+16] =           7307872110ull; c_[36*max_+17] =           8597496600ull;
      c_[36*max_+18] =           9075135300ull; c_[36*max_+19] =           8597496600ull; c_[36*max_+20] =           7307872110ull; c_[36*max_+21] =           5567902560ull;
      c_[36*max_+22] =           3796297200ull; c_[36*max_+23] =           2310789600ull; c_[36*max_+24] =           1251677700ull; c_[36*max_+25] =            600805296ull;
      c_[36*max_+26] =            254186856ull; c_[36*max_+27] =             94143280ull; c_[36*max_+28] =             30260340ull; c_[36*max_+29] =              8347680ull;
      c_[36*max_+30] =              1947792ull; c_[36*max_+31] =               376992ull; c_[36*max_+32] =                58905ull; c_[36*max_+33] =                 7140ull;
      c_[36*max_+34] =                  630ull; c_[36*max_+35] =                   36ull; c_[36*max_+36] =                    1ull; c_[37*max_+ 0] =                    1ull;
      c_[37*max_+ 1] =                   37ull; c_[37*max_+ 2] =                  666ull; c_[37*max_+ 3] =                 7770ull; c_[37*max_+ 4] =                66045ull;
      c_[37*max_+ 5] =               435897ull; c_[37*max_+ 6] =              2324784ull; c_[37*max_+ 7] =             10295472ull; c_[37*max_+ 8] =             38608020ull;
      c_[37*max_+ 9] =            124403620ull; c_[37*max_+10] =            348330136ull; c_[37*max_+11] =            854992152ull; c_[37*max_+12] =           1852482996ull;
      c_[37*max_+13] =           3562467300ull; c_[37*max_+14] =           6107086800ull; c_[37*max_+15] =           9364199760ull; c_[37*max_+16] =          12875774670ull;
      c_[37*max_+17] =          15905368710ull; c_[37*max_+18] =          17672631900ull; c_[37*max_+19] =          17672631900ull; c_[37*max_+20] =          15905368710ull;
      c_[37*max_+21] =          12875774670ull; c_[37*max_+22] =           9364199760ull; c_[37*max_+23] =           6107086800ull; c_[37*max_+24] =           3562467300ull;
      c_[37*max_+25] =           1852482996ull; c_[37*max_+26] =            854992152ull; c_[37*max_+27] =            348330136ull; c_[37*max_+28] =            124403620ull;
      c_[37*max_+29] =             38608020ull; c_[37*max_+30] =             10295472ull; c_[37*max_+31] =              2324784ull; c_[37*max_+32] =               435897ull;
      c_[37*max_+33] =                66045ull; c_[37*max_+34] =                 7770ull; c_[37*max_+35] =                  666ull; c_[37*max_+36] =                   37ull;
      c_[37*max_+37] =                    1ull; c_[38*max_+ 0] =                    1ull; c_[38*max_+ 1] =                   38ull; c_[38*max_+ 2] =                  703ull;
      c_[38*max_+ 3] =                 8436ull; c_[38*max_+ 4] =                73815ull; c_[38*max_+ 5] =               501942ull; c_[38*max_+ 6] =              2760681ull;
      c_[38*max_+ 7] =             12620256ull; c_[38*max_+ 8] =             48903492ull; c_[38*max_+ 9] =            163011640ull; c_[38*max_+10] =            472733756ull;
      c_[38*max_+11] =           1203322288ull; c_[38*max_+12] =           2707475148ull; c_[38*max_+13] =           5414950296ull; c_[38*max_+14] =           9669554100ull;
      c_[38*max_+15] =          15471286560ull; c_[38*max_+16] =          22239974430ull; c_[38*max_+17] =          28781143380ull; c_[38*max_+18] =          33578000610ull;
      c_[38*max_+19] =          35345263800ull; c_[38*max_+20] =          33578000610ull; c_[38*max_+21] =          28781143380ull; c_[38*max_+22] =          22239974430ull;
      c_[38*max_+23] =          15471286560ull; c_[38*max_+24] =           9669554100ull; c_[38*max_+25] =           5414950296ull; c_[38*max_+26] =           2707475148ull;
      c_[38*max_+27] =           1203322288ull; c_[38*max_+28] =            472733756ull; c_[38*max_+29] =            163011640ull; c_[38*max_+30] =             48903492ull;
      c_[38*max_+31] =             12620256ull; c_[38*max_+32] =              2760681ull; c_[38*max_+33] =               501942ull; c_[38*max_+34] =                73815ull;
      c_[38*max_+35] =                 8436ull; c_[38*max_+36] =                  703ull; c_[38*max_+37] =                   38ull; c_[38*max_+38] =                    1ull;
      c_[39*max_+ 0] =                    1ull; c_[39*max_+ 1] =                   39ull; c_[39*max_+ 2] =                  741ull; c_[39*max_+ 3] =                 9139ull;
      c_[39*max_+ 4] =                82251ull; c_[39*max_+ 5] =               575757ull; c_[39*max_+ 6] =              3262623ull; c_[39*max_+ 7] =             15380937ull;
      c_[39*max_+ 8] =             61523748ull; c_[39*max_+ 9] =            211915132ull; c_[39*max_+10] =            635745396ull; c_[39*max_+11] =           1676056044ull;
      c_[39*max_+12] =           3910797436ull; c_[39*max_+13] =           8122425444ull; c_[39*max_+14] =          15084504396ull; c_[39*max_+15] =          25140840660ull;
      c_[39*max_+16] =          37711260990ull; c_[39*max_+17] =          51021117810ull; c_[39*max_+18] =          62359143990ull; c_[39*max_+19] =          68923264410ull;
      c_[39*max_+20] =          68923264410ull; c_[39*max_+21] =          62359143990ull; c_[39*max_+22] =          51021117810ull; c_[39*max_+23] =          37711260990ull;
      c_[39*max_+24] =          25140840660ull; c_[39*max_+25] =          15084504396ull; c_[39*max_+26] =           8122425444ull; c_[39*max_+27] =           3910797436ull;
      c_[39*max_+28] =           1676056044ull; c_[39*max_+29] =            635745396ull; c_[39*max_+30] =            211915132ull; c_[39*max_+31] =             61523748ull;
      c_[39*max_+32] =             15380937ull; c_[39*max_+33] =              3262623ull; c_[39*max_+34] =               575757ull; c_[39*max_+35] =                82251ull;
      c_[39*max_+36] =                 9139ull; c_[39*max_+37] =                  741ull; c_[39*max_+38] =                   39ull; c_[39*max_+39] =                    1ull;
      c_[40*max_+ 0] =                    1ull; c_[40*max_+ 1] =                   40ull; c_[40*max_+ 2] =                  780ull; c_[40*max_+ 3] =                 9880ull;
      c_[40*max_+ 4] =                91390ull; c_[40*max_+ 5] =               658008ull; c_[40*max_+ 6] =              3838380ull; c_[40*max_+ 7] =             18643560ull;
      c_[40*max_+ 8] =             76904685ull; c_[40*max_+ 9] =            273438880ull; c_[40*max_+10] =            847660528ull; c_[40*max_+11] =           2311801440ull;
      c_[40*max_+12] =           5586853480ull; c_[40*max_+13] =          12033222880ull; c_[40*max_+14] =          23206929840ull; c_[40*max_+15] =          40225345056ull;
      c_[40*max_+16] =          62852101650ull; c_[40*max_+17] =          88732378800ull; c_[40*max_+18] =         113380261800ull; c_[40*max_+19] =         131282408400ull;
      c_[40*max_+20] =         137846528820ull; c_[40*max_+21] =         131282408400ull; c_[40*max_+22] =         113380261800ull; c_[40*max_+23] =          88732378800ull;
      c_[40*max_+24] =          62852101650ull; c_[40*max_+25] =          40225345056ull; c_[40*max_+26] =          23206929840ull; c_[40*max_+27] =          12033222880ull;
      c_[40*max_+28] =           5586853480ull; c_[40*max_+29] =           2311801440ull; c_[40*max_+30] =            847660528ull; c_[40*max_+31] =            273438880ull;
      c_[40*max_+32] =             76904685ull; c_[40*max_+33] =             18643560ull; c_[40*max_+34] =              3838380ull; c_[40*max_+35] =               658008ull;
      c_[40*max_+36] =                91390ull; c_[40*max_+37] =                 9880ull; c_[40*max_+38] =                  780ull; c_[40*max_+39] =                   40ull;
      c_[40*max_+40] =                    1ull; c_[41*max_+ 0] =                    1ull; c_[41*max_+ 1] =                   41ull; c_[41*max_+ 2] =                  820ull;
      c_[41*max_+ 3] =                10660ull; c_[41*max_+ 4] =               101270ull; c_[41*max_+ 5] =               749398ull; c_[41*max_+ 6] =              4496388ull;
      c_[41*max_+ 7] =             22481940ull; c_[41*max_+ 8] =             95548245ull; c_[41*max_+ 9] =            350343565ull; c_[41*max_+10] =           1121099408ull;
      c_[41*max_+11] =           3159461968ull; c_[41*max_+12] =           7898654920ull; c_[41*max_+13] =          17620076360ull; c_[41*max_+14] =          35240152720ull;
      c_[41*max_+15] =          63432274896ull; c_[41*max_+16] =         103077446706ull; c_[41*max_+17] =         151584480450ull; c_[41*max_+18] =         202112640600ull;
      c_[41*max_+19] =         244662670200ull; c_[41*max_+20] =         269128937220ull; c_[41*max_+21] =         269128937220ull; c_[41*max_+22] =         244662670200ull;
      c_[41*max_+23] =         202112640600ull; c_[41*max_+24] =         151584480450ull; c_[41*max_+25] =         103077446706ull; c_[41*max_+26] =          63432274896ull;
      c_[41*max_+27] =          35240152720ull; c_[41*max_+28] =          17620076360ull; c_[41*max_+29] =           7898654920ull; c_[41*max_+30] =           3159461968ull;
      c_[41*max_+31] =           1121099408ull; c_[41*max_+32] =            350343565ull; c_[41*max_+33] =             95548245ull; c_[41*max_+34] =             22481940ull;
      c_[41*max_+35] =              4496388ull; c_[41*max_+36] =               749398ull; c_[41*max_+37] =               101270ull; c_[41*max_+38] =                10660ull;
      c_[41*max_+39] =                  820ull; c_[41*max_+40] =                   41ull; c_[41*max_+41] =                    1ull; c_[42*max_+ 0] =                    1ull;
      c_[42*max_+ 1] =                   42ull; c_[42*max_+ 2] =                  861ull; c_[42*max_+ 3] =                11480ull; c_[42*max_+ 4] =               111930ull;
      c_[42*max_+ 5] =               850668ull; c_[42*max_+ 6] =              5245786ull; c_[42*max_+ 7] =             26978328ull; c_[42*max_+ 8] =            118030185ull;
      c_[42*max_+ 9] =            445891810ull; c_[42*max_+10] =           1471442973ull; c_[42*max_+11] =           4280561376ull; c_[42*max_+12] =          11058116888ull;
      c_[42*max_+13] =          25518731280ull; c_[42*max_+14] =          52860229080ull; c_[42*max_+15] =          98672427616ull; c_[42*max_+16] =         166509721602ull;
      c_[42*max_+17] =         254661927156ull; c_[42*max_+18] =         353697121050ull; c_[42*max_+19] =         446775310800ull; c_[42*max_+20] =         513791607420ull;
      c_[42*max_+21] =         538257874440ull; c_[42*max_+22] =         513791607420ull; c_[42*max_+23] =         446775310800ull; c_[42*max_+24] =         353697121050ull;
      c_[42*max_+25] =         254661927156ull; c_[42*max_+26] =         166509721602ull; c_[42*max_+27] =          98672427616ull; c_[42*max_+28] =          52860229080ull;
      c_[42*max_+29] =          25518731280ull; c_[42*max_+30] =          11058116888ull; c_[42*max_+31] =           4280561376ull; c_[42*max_+32] =           1471442973ull;
      c_[42*max_+33] =            445891810ull; c_[42*max_+34] =            118030185ull; c_[42*max_+35] =             26978328ull; c_[42*max_+36] =              5245786ull;
      c_[42*max_+37] =               850668ull; c_[42*max_+38] =               111930ull; c_[42*max_+39] =                11480ull; c_[42*max_+40] =                  861ull;
      c_[42*max_+41] =                   42ull; c_[42*max_+42] =                    1ull; c_[43*max_+ 0] =                    1ull; c_[43*max_+ 1] =                   43ull;
      c_[43*max_+ 2] =                  903ull; c_[43*max_+ 3] =                12341ull; c_[43*max_+ 4] =               123410ull; c_[43*max_+ 5] =               962598ull;
      c_[43*max_+ 6] =              6096454ull; c_[43*max_+ 7] =             32224114ull; c_[43*max_+ 8] =            145008513ull; c_[43*max_+ 9] =            563921995ull;
      c_[43*max_+10] =           1917334783ull; c_[43*max_+11] =           5752004349ull; c_[43*max_+12] =          15338678264ull; c_[43*max_+13] =          36576848168ull;
      c_[43*max_+14] =          78378960360ull; c_[43*max_+15] =         151532656696ull; c_[43*max_+16] =         265182149218ull; c_[43*max_+17] =         421171648758ull;
      c_[43*max_+18] =         608359048206ull; c_[43*max_+19] =         800472431850ull; c_[43*max_+20] =         960566918220ull; c_[43*max_+21] =        1052049481860ull;
      c_[43*max_+22] =        1052049481860ull; c_[43*max_+23] =         960566918220ull; c_[43*max_+24] =         800472431850ull; c_[43*max_+25] =         608359048206ull;
      c_[43*max_+26] =         421171648758ull; c_[43*max_+27] =         265182149218ull; c_[43*max_+28] =         151532656696ull; c_[43*max_+29] =          78378960360ull;
      c_[43*max_+30] =          36576848168ull; c_[43*max_+31] =          15338678264ull; c_[43*max_+32] =           5752004349ull; c_[43*max_+33] =           1917334783ull;
      c_[43*max_+34] =            563921995ull; c_[43*max_+35] =            145008513ull; c_[43*max_+36] =             32224114ull; c_[43*max_+37] =              6096454ull;
      c_[43*max_+38] =               962598ull; c_[43*max_+39] =               123410ull; c_[43*max_+40] =                12341ull; c_[43*max_+41] =                  903ull;
      c_[43*max_+42] =                   43ull; c_[43*max_+43] =                    1ull; c_[44*max_+ 0] =                    1ull; c_[44*max_+ 1] =                   44ull;
      c_[44*max_+ 2] =                  946ull; c_[44*max_+ 3] =                13244ull; c_[44*max_+ 4] =               135751ull; c_[44*max_+ 5] =              1086008ull;
      c_[44*max_+ 6] =              7059052ull; c_[44*max_+ 7] =             38320568ull; c_[44*max_+ 8] =            177232627ull; c_[44*max_+ 9] =            708930508ull;
      c_[44*max_+10] =           2481256778ull; c_[44*max_+11] =           7669339132ull; c_[44*max_+12] =          21090682613ull; c_[44*max_+13] =          51915526432ull;
      c_[44*max_+14] =         114955808528ull; c_[44*max_+15] =         229911617056ull; c_[44*max_+16] =         416714805914ull; c_[44*max_+17] =         686353797976ull;
      c_[44*max_+18] =        1029530696964ull; c_[44*max_+19] =        1408831480056ull; c_[44*max_+20] =        1761039350070ull; c_[44*max_+21] =        2012616400080ull;
      c_[44*max_+22] =        2104098963720ull; c_[44*max_+23] =        2012616400080ull; c_[44*max_+24] =        1761039350070ull; c_[44*max_+25] =        1408831480056ull;
      c_[44*max_+26] =        1029530696964ull; c_[44*max_+27] =         686353797976ull; c_[44*max_+28] =         416714805914ull; c_[44*max_+29] =         229911617056ull;
      c_[44*max_+30] =         114955808528ull; c_[44*max_+31] =          51915526432ull; c_[44*max_+32] =          21090682613ull; c_[44*max_+33] =           7669339132ull;
      c_[44*max_+34] =           2481256778ull; c_[44*max_+35] =            708930508ull; c_[44*max_+36] =            177232627ull; c_[44*max_+37] =             38320568ull;
      c_[44*max_+38] =              7059052ull; c_[44*max_+39] =              1086008ull; c_[44*max_+40] =               135751ull; c_[44*max_+41] =                13244ull;
      c_[44*max_+42] =                  946ull; c_[44*max_+43] =                   44ull; c_[44*max_+44] =                    1ull; c_[45*max_+ 0] =                    1ull;
      c_[45*max_+ 1] =                   45ull; c_[45*max_+ 2] =                  990ull; c_[45*max_+ 3] =                14190ull; c_[45*max_+ 4] =               148995ull;
      c_[45*max_+ 5] =              1221759ull; c_[45*max_+ 6] =              8145060ull; c_[45*max_+ 7] =             45379620ull; c_[45*max_+ 8] =            215553195ull;
      c_[45*max_+ 9] =            886163135ull; c_[45*max_+10] =           3190187286ull; c_[45*max_+11] =          10150595910ull; c_[45*max_+12] =          28760021745ull;
      c_[45*max_+13] =          73006209045ull; c_[45*max_+14] =         166871334960ull; c_[45*max_+15] =         344867425584ull; c_[45*max_+16] =         646626422970ull;
      c_[45*max_+17] =        1103068603890ull; c_[45*max_+18] =        1715884494940ull; c_[45*max_+19] =        2438362177020ull; c_[45*max_+20] =        3169870830126ull;
      c_[45*max_+21] =        3773655750150ull; c_[45*max_+22] =        4116715363800ull; c_[45*max_+23] =        4116715363800ull; c_[45*max_+24] =        3773655750150ull;
      c_[45*max_+25] =        3169870830126ull; c_[45*max_+26] =        2438362177020ull; c_[45*max_+27] =        1715884494940ull; c_[45*max_+28] =        1103068603890ull;
      c_[45*max_+29] =         646626422970ull; c_[45*max_+30] =         344867425584ull; c_[45*max_+31] =         166871334960ull; c_[45*max_+32] =          73006209045ull;
      c_[45*max_+33] =          28760021745ull; c_[45*max_+34] =          10150595910ull; c_[45*max_+35] =           3190187286ull; c_[45*max_+36] =            886163135ull;
      c_[45*max_+37] =            215553195ull; c_[45*max_+38] =             45379620ull; c_[45*max_+39] =              8145060ull; c_[45*max_+40] =              1221759ull;
      c_[45*max_+41] =               148995ull; c_[45*max_+42] =                14190ull; c_[45*max_+43] =                  990ull; c_[45*max_+44] =                   45ull;
      c_[45*max_+45] =                    1ull; c_[46*max_+ 0] =                    1ull; c_[46*max_+ 1] =                   46ull; c_[46*max_+ 2] =                 1035ull;
      c_[46*max_+ 3] =                15180ull; c_[46*max_+ 4] =               163185ull; c_[46*max_+ 5] =              1370754ull; c_[46*max_+ 6] =              9366819ull;
      c_[46*max_+ 7] =             53524680ull; c_[46*max_+ 8] =            260932815ull; c_[46*max_+ 9] =           1101716330ull; c_[46*max_+10] =           4076350421ull;
      c_[46*max_+11] =          13340783196ull; c_[46*max_+12] =          38910617655ull; c_[46*max_+13] =         101766230790ull; c_[46*max_+14] =         239877544005ull;
      c_[46*max_+15] =         511738760544ull; c_[46*max_+16] =         991493848554ull; c_[46*max_+17] =        1749695026860ull; c_[46*max_+18] =        2818953098830ull;
      c_[46*max_+19] =        4154246671960ull; c_[46*max_+20] =        5608233007146ull; c_[46*max_+21] =        6943526580276ull; c_[46*max_+22] =        7890371113950ull;
      c_[46*max_+23] =        8233430727600ull; c_[46*max_+24] =        7890371113950ull; c_[46*max_+25] =        6943526580276ull; c_[46*max_+26] =        5608233007146ull;
      c_[46*max_+27] =        4154246671960ull; c_[46*max_+28] =        2818953098830ull; c_[46*max_+29] =        1749695026860ull; c_[46*max_+30] =         991493848554ull;
      c_[46*max_+31] =         511738760544ull; c_[46*max_+32] =         239877544005ull; c_[46*max_+33] =         101766230790ull; c_[46*max_+34] =          38910617655ull;
      c_[46*max_+35] =          13340783196ull; c_[46*max_+36] =           4076350421ull; c_[46*max_+37] =           1101716330ull; c_[46*max_+38] =            260932815ull;
      c_[46*max_+39] =             53524680ull; c_[46*max_+40] =              9366819ull; c_[46*max_+41] =              1370754ull; c_[46*max_+42] =               163185ull;
      c_[46*max_+43] =                15180ull; c_[46*max_+44] =                 1035ull; c_[46*max_+45] =                   46ull; c_[46*max_+46] =                    1ull;
      c_[47*max_+ 0] =                    1ull; c_[47*max_+ 1] =                   47ull; c_[47*max_+ 2] =                 1081ull; c_[47*max_+ 3] =                16215ull;
      c_[47*max_+ 4] =               178365ull; c_[47*max_+ 5] =              1533939ull; c_[47*max_+ 6] =             10737573ull; c_[47*max_+ 7] =             62891499ull;
      c_[47*max_+ 8] =            314457495ull; c_[47*max_+ 9] =           1362649145ull; c_[47*max_+10] =           5178066751ull; c_[47*max_+11] =          17417133617ull;
      c_[47*max_+12] =          52251400851ull; c_[47*max_+13] =         140676848445ull; c_[47*max_+14] =         341643774795ull; c_[47*max_+15] =         751616304549ull;
      c_[47*max_+16] =        1503232609098ull; c_[47*max_+17] =        2741188875414ull; c_[47*max_+18] =        4568648125690ull; c_[47*max_+19] =        6973199770790ull;
      c_[47*max_+20] =        9762479679106ull; c_[47*max_+21] =       12551759587422ull; c_[47*max_+22] =       14833897694226ull; c_[47*max_+23] =       16123801841550ull;
      c_[47*max_+24] =       16123801841550ull; c_[47*max_+25] =       14833897694226ull; c_[47*max_+26] =       12551759587422ull; c_[47*max_+27] =        9762479679106ull;
      c_[47*max_+28] =        6973199770790ull; c_[47*max_+29] =        4568648125690ull; c_[47*max_+30] =        2741188875414ull; c_[47*max_+31] =        1503232609098ull;
      c_[47*max_+32] =         751616304549ull; c_[47*max_+33] =         341643774795ull; c_[47*max_+34] =         140676848445ull; c_[47*max_+35] =          52251400851ull;
      c_[47*max_+36] =          17417133617ull; c_[47*max_+37] =           5178066751ull; c_[47*max_+38] =           1362649145ull; c_[47*max_+39] =            314457495ull;
      c_[47*max_+40] =             62891499ull; c_[47*max_+41] =             10737573ull; c_[47*max_+42] =              1533939ull; c_[47*max_+43] =               178365ull;
      c_[47*max_+44] =                16215ull; c_[47*max_+45] =                 1081ull; c_[47*max_+46] =                   47ull; c_[47*max_+47] =                    1ull;
      c_[48*max_+ 0] =                    1ull; c_[48*max_+ 1] =                   48ull; c_[48*max_+ 2] =                 1128ull; c_[48*max_+ 3] =                17296ull;
      c_[48*max_+ 4] =               194580ull; c_[48*max_+ 5] =              1712304ull; c_[48*max_+ 6] =             12271512ull; c_[48*max_+ 7] =             73629072ull;
      c_[48*max_+ 8] =            377348994ull; c_[48*max_+ 9] =           1677106640ull; c_[48*max_+10] =           6540715896ull; c_[48*max_+11] =          22595200368ull;
      c_[48*max_+12] =          69668534468ull; c_[48*max_+13] =         192928249296ull; c_[48*max_+14] =         482320623240ull; c_[48*max_+15] =        1093260079344ull;
      c_[48*max_+16] =        2254848913647ull; c_[48*max_+17] =        4244421484512ull; c_[48*max_+18] =        7309837001104ull; c_[48*max_+19] =       11541847896480ull;
      c_[48*max_+20] =       16735679449896ull; c_[48*max_+21] =       22314239266528ull; c_[48*max_+22] =       27385657281648ull; c_[48*max_+23] =       30957699535776ull;
      c_[48*max_+24] =       32247603683100ull; c_[48*max_+25] =       30957699535776ull; c_[48*max_+26] =       27385657281648ull; c_[48*max_+27] =       22314239266528ull;
      c_[48*max_+28] =       16735679449896ull; c_[48*max_+29] =       11541847896480ull; c_[48*max_+30] =        7309837001104ull; c_[48*max_+31] =        4244421484512ull;
      c_[48*max_+32] =        2254848913647ull; c_[48*max_+33] =        1093260079344ull; c_[48*max_+34] =         482320623240ull; c_[48*max_+35] =         192928249296ull;
      c_[48*max_+36] =          69668534468ull; c_[48*max_+37] =          22595200368ull; c_[48*max_+38] =           6540715896ull; c_[48*max_+39] =           1677106640ull;
      c_[48*max_+40] =            377348994ull; c_[48*max_+41] =             73629072ull; c_[48*max_+42] =             12271512ull; c_[48*max_+43] =              1712304ull;
      c_[48*max_+44] =               194580ull; c_[48*max_+45] =                17296ull; c_[48*max_+46] =                 1128ull; c_[48*max_+47] =                   48ull;
      c_[48*max_+48] =                    1ull; c_[49*max_+ 0] =                    1ull; c_[49*max_+ 1] =                   49ull; c_[49*max_+ 2] =                 1176ull;
      c_[49*max_+ 3] =                18424ull; c_[49*max_+ 4] =               211876ull; c_[49*max_+ 5] =              1906884ull; c_[49*max_+ 6] =             13983816ull;
      c_[49*max_+ 7] =             85900584ull; c_[49*max_+ 8] =            450978066ull; c_[49*max_+ 9] =           2054455634ull; c_[49*max_+10] =           8217822536ull;
      c_[49*max_+11] =          29135916264ull; c_[49*max_+12] =          92263734836ull; c_[49*max_+13] =         262596783764ull; c_[49*max_+14] =         675248872536ull;
      c_[49*max_+15] =        1575580702584ull; c_[49*max_+16] =        3348108992991ull; c_[49*max_+17] =        6499270398159ull; c_[49*max_+18] =       11554258485616ull;
      c_[49*max_+19] =       18851684897584ull; c_[49*max_+20] =       28277527346376ull; c_[49*max_+21] =       39049918716424ull; c_[49*max_+22] =       49699896548176ull;
      c_[49*max_+23] =       58343356817424ull; c_[49*max_+24] =       63205303218876ull; c_[49*max_+25] =       63205303218876ull; c_[49*max_+26] =       58343356817424ull;
      c_[49*max_+27] =       49699896548176ull; c_[49*max_+28] =       39049918716424ull; c_[49*max_+29] =       28277527346376ull; c_[49*max_+30] =       18851684897584ull;
      c_[49*max_+31] =       11554258485616ull; c_[49*max_+32] =        6499270398159ull; c_[49*max_+33] =        3348108992991ull; c_[49*max_+34] =        1575580702584ull;
      c_[49*max_+35] =         675248872536ull; c_[49*max_+36] =         262596783764ull; c_[49*max_+37] =          92263734836ull; c_[49*max_+38] =          29135916264ull;
      c_[49*max_+39] =           8217822536ull; c_[49*max_+40] =           2054455634ull; c_[49*max_+41] =            450978066ull; c_[49*max_+42] =             85900584ull;
      c_[49*max_+43] =             13983816ull; c_[49*max_+44] =              1906884ull; c_[49*max_+45] =               211876ull; c_[49*max_+46] =                18424ull;
      c_[49*max_+47] =                 1176ull; c_[49*max_+48] =                   49ull; c_[49*max_+49] =                    1ull; c_[50*max_+ 0] =                    1ull;
      c_[50*max_+ 1] =                   50ull; c_[50*max_+ 2] =                 1225ull; c_[50*max_+ 3] =                19600ull; c_[50*max_+ 4] =               230300ull;
      c_[50*max_+ 5] =              2118760ull; c_[50*max_+ 6] =             15890700ull; c_[50*max_+ 7] =             99884400ull; c_[50*max_+ 8] =            536878650ull;
      c_[50*max_+ 9] =           2505433700ull; c_[50*max_+10] =          10272278170ull; c_[50*max_+11] =          37353738800ull; c_[50*max_+12] =         121399651100ull;
      c_[50*max_+13] =         354860518600ull; c_[50*max_+14] =         937845656300ull; c_[50*max_+15] =        2250829575120ull; c_[50*max_+16] =        4923689695575ull;
      c_[50*max_+17] =        9847379391150ull; c_[50*max_+18] =       18053528883775ull; c_[50*max_+19] =       30405943383200ull; c_[50*max_+20] =       47129212243960ull;
      c_[50*max_+21] =       67327446062800ull; c_[50*max_+22] =       88749815264600ull; c_[50*max_+23] =      108043253365600ull; c_[50*max_+24] =      121548660036300ull;
      c_[50*max_+25] =      126410606437752ull; c_[50*max_+26] =      121548660036300ull; c_[50*max_+27] =      108043253365600ull; c_[50*max_+28] =       88749815264600ull;
      c_[50*max_+29] =       67327446062800ull; c_[50*max_+30] =       47129212243960ull; c_[50*max_+31] =       30405943383200ull; c_[50*max_+32] =       18053528883775ull;
      c_[50*max_+33] =        9847379391150ull; c_[50*max_+34] =        4923689695575ull; c_[50*max_+35] =        2250829575120ull; c_[50*max_+36] =         937845656300ull;
      c_[50*max_+37] =         354860518600ull; c_[50*max_+38] =         121399651100ull; c_[50*max_+39] =          37353738800ull; c_[50*max_+40] =          10272278170ull;
      c_[50*max_+41] =           2505433700ull; c_[50*max_+42] =            536878650ull; c_[50*max_+43] =             99884400ull; c_[50*max_+44] =             15890700ull;
      c_[50*max_+45] =              2118760ull; c_[50*max_+46] =               230300ull; c_[50*max_+47] =                19600ull; c_[50*max_+48] =                 1225ull;
      c_[50*max_+49] =                   50ull; c_[50*max_+50] =                    1ull; c_[51*max_+ 0] =                    1ull; c_[51*max_+ 1] =                   51ull;
      c_[51*max_+ 2] =                 1275ull; c_[51*max_+ 3] =                20825ull; c_[51*max_+ 4] =               249900ull; c_[51*max_+ 5] =              2349060ull;
      c_[51*max_+ 6] =             18009460ull; c_[51*max_+ 7] =            115775100ull; c_[51*max_+ 8] =            636763050ull; c_[51*max_+ 9] =           3042312350ull;
      c_[51*max_+10] =          12777711870ull; c_[51*max_+11] =          47626016970ull; c_[51*max_+12] =         158753389900ull; c_[51*max_+13] =         476260169700ull;
      c_[51*max_+14] =        1292706174900ull; c_[51*max_+15] =        3188675231420ull; c_[51*max_+16] =        7174519270695ull; c_[51*max_+17] =       14771069086725ull;
      c_[51*max_+18] =       27900908274925ull; c_[51*max_+19] =       48459472266975ull; c_[51*max_+20] =       77535155627160ull; c_[51*max_+21] =      114456658306760ull;
      c_[51*max_+22] =      156077261327400ull; c_[51*max_+23] =      196793068630200ull; c_[51*max_+24] =      229591913401900ull; c_[51*max_+25] =      247959266474052ull;
      c_[51*max_+26] =      247959266474052ull; c_[51*max_+27] =      229591913401900ull; c_[51*max_+28] =      196793068630200ull; c_[51*max_+29] =      156077261327400ull;
      c_[51*max_+30] =      114456658306760ull; c_[51*max_+31] =       77535155627160ull; c_[51*max_+32] =       48459472266975ull; c_[51*max_+33] =       27900908274925ull;
      c_[51*max_+34] =       14771069086725ull; c_[51*max_+35] =        7174519270695ull; c_[51*max_+36] =        3188675231420ull; c_[51*max_+37] =        1292706174900ull;
      c_[51*max_+38] =         476260169700ull; c_[51*max_+39] =         158753389900ull; c_[51*max_+40] =          47626016970ull; c_[51*max_+41] =          12777711870ull;
      c_[51*max_+42] =           3042312350ull; c_[51*max_+43] =            636763050ull; c_[51*max_+44] =            115775100ull; c_[51*max_+45] =             18009460ull;
      c_[51*max_+46] =              2349060ull; c_[51*max_+47] =               249900ull; c_[51*max_+48] =                20825ull; c_[51*max_+49] =                 1275ull;
      c_[51*max_+50] =                   51ull; c_[51*max_+51] =                    1ull; c_[52*max_+ 0] =                    1ull; c_[52*max_+ 1] =                   52ull;
      c_[52*max_+ 2] =                 1326ull; c_[52*max_+ 3] =                22100ull; c_[52*max_+ 4] =               270725ull; c_[52*max_+ 5] =              2598960ull;
      c_[52*max_+ 6] =             20358520ull; c_[52*max_+ 7] =            133784560ull; c_[52*max_+ 8] =            752538150ull; c_[52*max_+ 9] =           3679075400ull;
      c_[52*max_+10] =          15820024220ull; c_[52*max_+11] =          60403728840ull; c_[52*max_+12] =         206379406870ull; c_[52*max_+13] =         635013559600ull;
      c_[52*max_+14] =        1768966344600ull; c_[52*max_+15] =        4481381406320ull; c_[52*max_+16] =       10363194502115ull; c_[52*max_+17] =       21945588357420ull;
      c_[52*max_+18] =       42671977361650ull; c_[52*max_+19] =       76360380541900ull; c_[52*max_+20] =      125994627894135ull; c_[52*max_+21] =      191991813933920ull;
      c_[52*max_+22] =      270533919634160ull; c_[52*max_+23] =      352870329957600ull; c_[52*max_+24] =      426384982032100ull; c_[52*max_+25] =      477551179875952ull;
      c_[52*max_+26] =      495918532948104ull; c_[52*max_+27] =      477551179875952ull; c_[52*max_+28] =      426384982032100ull; c_[52*max_+29] =      352870329957600ull;
      c_[52*max_+30] =      270533919634160ull; c_[52*max_+31] =      191991813933920ull; c_[52*max_+32] =      125994627894135ull; c_[52*max_+33] =       76360380541900ull;
      c_[52*max_+34] =       42671977361650ull; c_[52*max_+35] =       21945588357420ull; c_[52*max_+36] =       10363194502115ull; c_[52*max_+37] =        4481381406320ull;
      c_[52*max_+38] =        1768966344600ull; c_[52*max_+39] =         635013559600ull; c_[52*max_+40] =         206379406870ull; c_[52*max_+41] =          60403728840ull;
      c_[52*max_+42] =          15820024220ull; c_[52*max_+43] =           3679075400ull; c_[52*max_+44] =            752538150ull; c_[52*max_+45] =            133784560ull;
      c_[52*max_+46] =             20358520ull; c_[52*max_+47] =              2598960ull; c_[52*max_+48] =               270725ull; c_[52*max_+49] =                22100ull;
      c_[52*max_+50] =                 1326ull; c_[52*max_+51] =                   52ull; c_[52*max_+52] =                    1ull; c_[53*max_+ 0] =                    1ull;
      c_[53*max_+ 1] =                   53ull; c_[53*max_+ 2] =                 1378ull; c_[53*max_+ 3] =                23426ull; c_[53*max_+ 4] =               292825ull;
      c_[53*max_+ 5] =              2869685ull; c_[53*max_+ 6] =             22957480ull; c_[53*max_+ 7] =            154143080ull; c_[53*max_+ 8] =            886322710ull;
      c_[53*max_+ 9] =           4431613550ull; c_[53*max_+10] =          19499099620ull; c_[53*max_+11] =          76223753060ull; c_[53*max_+12] =         266783135710ull;
      c_[53*max_+13] =         841392966470ull; c_[53*max_+14] =        2403979904200ull; c_[53*max_+15] =        6250347750920ull; c_[53*max_+16] =       14844575908435ull;
      c_[53*max_+17] =       32308782859535ull; c_[53*max_+18] =       64617565719070ull; c_[53*max_+19] =      119032357903550ull; c_[53*max_+20] =      202355008436035ull;
      c_[53*max_+21] =      317986441828055ull; c_[53*max_+22] =      462525733568080ull; c_[53*max_+23] =      623404249591760ull; c_[53*max_+24] =      779255311989700ull;
      c_[53*max_+25] =      903936161908052ull; c_[53*max_+26] =      973469712824056ull; c_[53*max_+27] =      973469712824056ull; c_[53*max_+28] =      903936161908052ull;
      c_[53*max_+29] =      779255311989700ull; c_[53*max_+30] =      623404249591760ull; c_[53*max_+31] =      462525733568080ull; c_[53*max_+32] =      317986441828055ull;
      c_[53*max_+33] =      202355008436035ull; c_[53*max_+34] =      119032357903550ull; c_[53*max_+35] =       64617565719070ull; c_[53*max_+36] =       32308782859535ull;
      c_[53*max_+37] =       14844575908435ull; c_[53*max_+38] =        6250347750920ull; c_[53*max_+39] =        2403979904200ull; c_[53*max_+40] =         841392966470ull;
      c_[53*max_+41] =         266783135710ull; c_[53*max_+42] =          76223753060ull; c_[53*max_+43] =          19499099620ull; c_[53*max_+44] =           4431613550ull;
      c_[53*max_+45] =            886322710ull; c_[53*max_+46] =            154143080ull; c_[53*max_+47] =             22957480ull; c_[53*max_+48] =              2869685ull;
      c_[53*max_+49] =               292825ull; c_[53*max_+50] =                23426ull; c_[53*max_+51] =                 1378ull; c_[53*max_+52] =                   53ull;
      c_[53*max_+53] =                    1ull; c_[54*max_+ 0] =                    1ull; c_[54*max_+ 1] =                   54ull; c_[54*max_+ 2] =                 1431ull;
      c_[54*max_+ 3] =                24804ull; c_[54*max_+ 4] =               316251ull; c_[54*max_+ 5] =              3162510ull; c_[54*max_+ 6] =             25827165ull;
      c_[54*max_+ 7] =            177100560ull; c_[54*max_+ 8] =           1040465790ull; c_[54*max_+ 9] =           5317936260ull; c_[54*max_+10] =          23930713170ull;
      c_[54*max_+11] =          95722852680ull; c_[54*max_+12] =         343006888770ull; c_[54*max_+13] =        1108176102180ull; c_[54*max_+14] =        3245372870670ull;
      c_[54*max_+15] =        8654327655120ull; c_[54*max_+16] =       21094923659355ull; c_[54*max_+17] =       47153358767970ull; c_[54*max_+18] =       96926348578605ull;
      c_[54*max_+19] =      183649923622620ull; c_[54*max_+20] =      321387366339585ull; c_[54*max_+21] =      520341450264090ull; c_[54*max_+22] =      780512175396135ull;
      c_[54*max_+23] =     1085929983159840ull; c_[54*max_+24] =     1402659561581460ull; c_[54*max_+25] =     1683191473897752ull; c_[54*max_+26] =     1877405874732108ull;
      c_[54*max_+27] =     1946939425648112ull; c_[54*max_+28] =     1877405874732108ull; c_[54*max_+29] =     1683191473897752ull; c_[54*max_+30] =     1402659561581460ull;
      c_[54*max_+31] =     1085929983159840ull; c_[54*max_+32] =      780512175396135ull; c_[54*max_+33] =      520341450264090ull; c_[54*max_+34] =      321387366339585ull;
      c_[54*max_+35] =      183649923622620ull; c_[54*max_+36] =       96926348578605ull; c_[54*max_+37] =       47153358767970ull; c_[54*max_+38] =       21094923659355ull;
      c_[54*max_+39] =        8654327655120ull; c_[54*max_+40] =        3245372870670ull; c_[54*max_+41] =        1108176102180ull; c_[54*max_+42] =         343006888770ull;
      c_[54*max_+43] =          95722852680ull; c_[54*max_+44] =          23930713170ull; c_[54*max_+45] =           5317936260ull; c_[54*max_+46] =           1040465790ull;
      c_[54*max_+47] =            177100560ull; c_[54*max_+48] =             25827165ull; c_[54*max_+49] =              3162510ull; c_[54*max_+50] =               316251ull;
      c_[54*max_+51] =                24804ull; c_[54*max_+52] =                 1431ull; c_[54*max_+53] =                   54ull; c_[54*max_+54] =                    1ull;
      c_[55*max_+ 0] =                    1ull; c_[55*max_+ 1] =                   55ull; c_[55*max_+ 2] =                 1485ull; c_[55*max_+ 3] =                26235ull;
      c_[55*max_+ 4] =               341055ull; c_[55*max_+ 5] =              3478761ull; c_[55*max_+ 6] =             28989675ull; c_[55*max_+ 7] =            202927725ull;
      c_[55*max_+ 8] =           1217566350ull; c_[55*max_+ 9] =           6358402050ull; c_[55*max_+10] =          29248649430ull; c_[55*max_+11] =         119653565850ull;
      c_[55*max_+12] =         438729741450ull; c_[55*max_+13] =        1451182990950ull; c_[55*max_+14] =        4353548972850ull; c_[55*max_+15] =       11899700525790ull;
      c_[55*max_+16] =       29749251314475ull; c_[55*max_+17] =       68248282427325ull; c_[55*max_+18] =      144079707346575ull; c_[55*max_+19] =      280576272201225ull;
      c_[55*max_+20] =      505037289962205ull; c_[55*max_+21] =      841728816603675ull; c_[55*max_+22] =     1300853625660225ull; c_[55*max_+23] =     1866442158555975ull;
      c_[55*max_+24] =     2488589544741300ull; c_[55*max_+25] =     3085851035479212ull; c_[55*max_+26] =     3560597348629860ull; c_[55*max_+27] =     3824345300380220ull;
      c_[55*max_+28] =     3824345300380220ull; c_[55*max_+29] =     3560597348629860ull; c_[55*max_+30] =     3085851035479212ull; c_[55*max_+31] =     2488589544741300ull;
      c_[55*max_+32] =     1866442158555975ull; c_[55*max_+33] =     1300853625660225ull; c_[55*max_+34] =      841728816603675ull; c_[55*max_+35] =      505037289962205ull;
      c_[55*max_+36] =      280576272201225ull; c_[55*max_+37] =      144079707346575ull; c_[55*max_+38] =       68248282427325ull; c_[55*max_+39] =       29749251314475ull;
      c_[55*max_+40] =       11899700525790ull; c_[55*max_+41] =        4353548972850ull; c_[55*max_+42] =        1451182990950ull; c_[55*max_+43] =         438729741450ull;
      c_[55*max_+44] =         119653565850ull; c_[55*max_+45] =          29248649430ull; c_[55*max_+46] =           6358402050ull; c_[55*max_+47] =           1217566350ull;
      c_[55*max_+48] =            202927725ull; c_[55*max_+49] =             28989675ull; c_[55*max_+50] =              3478761ull; c_[55*max_+51] =               341055ull;
      c_[55*max_+52] =                26235ull; c_[55*max_+53] =                 1485ull; c_[55*max_+54] =                   55ull; c_[55*max_+55] =                    1ull;
      c_[56*max_+ 0] =                    1ull; c_[56*max_+ 1] =                   56ull; c_[56*max_+ 2] =                 1540ull; c_[56*max_+ 3] =                27720ull;
      c_[56*max_+ 4] =               367290ull; c_[56*max_+ 5] =              3819816ull; c_[56*max_+ 6] =             32468436ull; c_[56*max_+ 7] =            231917400ull;
      c_[56*max_+ 8] =           1420494075ull; c_[56*max_+ 9] =           7575968400ull; c_[56*max_+10] =          35607051480ull; c_[56*max_+11] =         148902215280ull;
      c_[56*max_+12] =         558383307300ull; c_[56*max_+13] =        1889912732400ull; c_[56*max_+14] =        5804731963800ull; c_[56*max_+15] =       16253249498640ull;
      c_[56*max_+16] =       41648951840265ull; c_[56*max_+17] =       97997533741800ull; c_[56*max_+18] =      212327989773900ull; c_[56*max_+19] =      424655979547800ull;
      c_[56*max_+20] =      785613562163430ull; c_[56*max_+21] =     1346766106565880ull; c_[56*max_+22] =     2142582442263900ull; c_[56*max_+23] =     3167295784216200ull;
      c_[56*max_+24] =     4355031703297275ull; c_[56*max_+25] =     5574440580220512ull; c_[56*max_+26] =     6646448384109072ull; c_[56*max_+27] =     7384942649010080ull;
      c_[56*max_+28] =     7648690600760440ull; c_[56*max_+29] =     7384942649010080ull; c_[56*max_+30] =     6646448384109072ull; c_[56*max_+31] =     5574440580220512ull;
      c_[56*max_+32] =     4355031703297275ull; c_[56*max_+33] =     3167295784216200ull; c_[56*max_+34] =     2142582442263900ull; c_[56*max_+35] =     1346766106565880ull;
      c_[56*max_+36] =      785613562163430ull; c_[56*max_+37] =      424655979547800ull; c_[56*max_+38] =      212327989773900ull; c_[56*max_+39] =       97997533741800ull;
      c_[56*max_+40] =       41648951840265ull; c_[56*max_+41] =       16253249498640ull; c_[56*max_+42] =        5804731963800ull; c_[56*max_+43] =        1889912732400ull;
      c_[56*max_+44] =         558383307300ull; c_[56*max_+45] =         148902215280ull; c_[56*max_+46] =          35607051480ull; c_[56*max_+47] =           7575968400ull;
      c_[56*max_+48] =           1420494075ull; c_[56*max_+49] =            231917400ull; c_[56*max_+50] =             32468436ull; c_[56*max_+51] =              3819816ull;
      c_[56*max_+52] =               367290ull; c_[56*max_+53] =                27720ull; c_[56*max_+54] =                 1540ull; c_[56*max_+55] =                   56ull;
      c_[56*max_+56] =                    1ull; c_[57*max_+ 0] =                    1ull; c_[57*max_+ 1] =                   57ull; c_[57*max_+ 2] =                 1596ull;
      c_[57*max_+ 3] =                29260ull; c_[57*max_+ 4] =               395010ull; c_[57*max_+ 5] =              4187106ull; c_[57*max_+ 6] =             36288252ull;
      c_[57*max_+ 7] =            264385836ull; c_[57*max_+ 8] =           1652411475ull; c_[57*max_+ 9] =           8996462475ull; c_[57*max_+10] =          43183019880ull;
      c_[57*max_+11] =         184509266760ull; c_[57*max_+12] =         707285522580ull; c_[57*max_+13] =        2448296039700ull; c_[57*max_+14] =        7694644696200ull;
      c_[57*max_+15] =       22057981462440ull; c_[57*max_+16] =       57902201338905ull; c_[57*max_+17] =      139646485582065ull; c_[57*max_+18] =      310325523515700ull;
      c_[57*max_+19] =      636983969321700ull; c_[57*max_+20] =     1210269541711230ull; c_[57*max_+21] =     2132379668729310ull; c_[57*max_+22] =     3489348548829780ull;
      c_[57*max_+23] =     5309878226480100ull; c_[57*max_+24] =     7522327487513475ull; c_[57*max_+25] =     9929472283517787ull; c_[57*max_+26] =    12220888964329584ull;
      c_[57*max_+27] =    14031391033119152ull; c_[57*max_+28] =    15033633249770520ull; c_[57*max_+29] =    15033633249770520ull; c_[57*max_+30] =    14031391033119152ull;
      c_[57*max_+31] =    12220888964329584ull; c_[57*max_+32] =     9929472283517787ull; c_[57*max_+33] =     7522327487513475ull; c_[57*max_+34] =     5309878226480100ull;
      c_[57*max_+35] =     3489348548829780ull; c_[57*max_+36] =     2132379668729310ull; c_[57*max_+37] =     1210269541711230ull; c_[57*max_+38] =      636983969321700ull;
      c_[57*max_+39] =      310325523515700ull; c_[57*max_+40] =      139646485582065ull; c_[57*max_+41] =       57902201338905ull; c_[57*max_+42] =       22057981462440ull;
      c_[57*max_+43] =        7694644696200ull; c_[57*max_+44] =        2448296039700ull; c_[57*max_+45] =         707285522580ull; c_[57*max_+46] =         184509266760ull;
      c_[57*max_+47] =          43183019880ull; c_[57*max_+48] =           8996462475ull; c_[57*max_+49] =           1652411475ull; c_[57*max_+50] =            264385836ull;
      c_[57*max_+51] =             36288252ull; c_[57*max_+52] =              4187106ull; c_[57*max_+53] =               395010ull; c_[57*max_+54] =                29260ull;
      c_[57*max_+55] =                 1596ull; c_[57*max_+56] =                   57ull; c_[57*max_+57] =                    1ull; c_[58*max_+ 0] =                    1ull;
      c_[58*max_+ 1] =                   58ull; c_[58*max_+ 2] =                 1653ull; c_[58*max_+ 3] =                30856ull; c_[58*max_+ 4] =               424270ull;
      c_[58*max_+ 5] =              4582116ull; c_[58*max_+ 6] =             40475358ull; c_[58*max_+ 7] =            300674088ull; c_[58*max_+ 8] =           1916797311ull;
      c_[58*max_+ 9] =          10648873950ull; c_[58*max_+10] =          52179482355ull; c_[58*max_+11] =         227692286640ull; c_[58*max_+12] =         891794789340ull;
      c_[58*max_+13] =        3155581562280ull; c_[58*max_+14] =       10142940735900ull; c_[58*max_+15] =       29752626158640ull; c_[58*max_+16] =       79960182801345ull;
      c_[58*max_+17] =      197548686920970ull; c_[58*max_+18] =      449972009097765ull; c_[58*max_+19] =      947309492837400ull; c_[58*max_+20] =     1847253511032930ull;
      c_[58*max_+21] =     3342649210440540ull; c_[58*max_+22] =     5621728217559090ull; c_[58*max_+23] =     8799226775309880ull; c_[58*max_+24] =    12832205713993575ull;
      c_[58*max_+25] =    17451799771031262ull; c_[58*max_+26] =    22150361247847371ull; c_[58*max_+27] =    26252279997448736ull; c_[58*max_+28] =    29065024282889672ull;
      c_[58*max_+29] =    30067266499541040ull; c_[58*max_+30] =    29065024282889672ull; c_[58*max_+31] =    26252279997448736ull; c_[58*max_+32] =    22150361247847371ull;
      c_[58*max_+33] =    17451799771031262ull; c_[58*max_+34] =    12832205713993575ull; c_[58*max_+35] =     8799226775309880ull; c_[58*max_+36] =     5621728217559090ull;
      c_[58*max_+37] =     3342649210440540ull; c_[58*max_+38] =     1847253511032930ull; c_[58*max_+39] =      947309492837400ull; c_[58*max_+40] =      449972009097765ull;
      c_[58*max_+41] =      197548686920970ull; c_[58*max_+42] =       79960182801345ull; c_[58*max_+43] =       29752626158640ull; c_[58*max_+44] =       10142940735900ull;
      c_[58*max_+45] =        3155581562280ull; c_[58*max_+46] =         891794789340ull; c_[58*max_+47] =         227692286640ull; c_[58*max_+48] =          52179482355ull;
      c_[58*max_+49] =          10648873950ull; c_[58*max_+50] =           1916797311ull; c_[58*max_+51] =            300674088ull; c_[58*max_+52] =             40475358ull;
      c_[58*max_+53] =              4582116ull; c_[58*max_+54] =               424270ull; c_[58*max_+55] =                30856ull; c_[58*max_+56] =                 1653ull;
      c_[58*max_+57] =                   58ull; c_[58*max_+58] =                    1ull; c_[59*max_+ 0] =                    1ull; c_[59*max_+ 1] =                   59ull;
      c_[59*max_+ 2] =                 1711ull; c_[59*max_+ 3] =                32509ull; c_[59*max_+ 4] =               455126ull; c_[59*max_+ 5] =              5006386ull;
      c_[59*max_+ 6] =             45057474ull; c_[59*max_+ 7] =            341149446ull; c_[59*max_+ 8] =           2217471399ull; c_[59*max_+ 9] =          12565671261ull;
      c_[59*max_+10] =          62828356305ull; c_[59*max_+11] =         279871768995ull; c_[59*max_+12] =        1119487075980ull; c_[59*max_+13] =        4047376351620ull;
      c_[59*max_+14] =       13298522298180ull; c_[59*max_+15] =       39895566894540ull; c_[59*max_+16] =      109712808959985ull; c_[59*max_+17] =      277508869722315ull;
      c_[59*max_+18] =      647520696018735ull; c_[59*max_+19] =     1397281501935165ull; c_[59*max_+20] =     2794563003870330ull; c_[59*max_+21] =     5189902721473470ull;
      c_[59*max_+22] =     8964377427999630ull; c_[59*max_+23] =    14420954992868970ull; c_[59*max_+24] =    21631432489303455ull; c_[59*max_+25] =    30284005485024837ull;
      c_[59*max_+26] =    39602161018878633ull; c_[59*max_+27] =    48402641245296107ull; c_[59*max_+28] =    55317304280338408ull; c_[59*max_+29] =    59132290782430712ull;
      c_[59*max_+30] =    59132290782430712ull; c_[59*max_+31] =    55317304280338408ull; c_[59*max_+32] =    48402641245296107ull; c_[59*max_+33] =    39602161018878633ull;
      c_[59*max_+34] =    30284005485024837ull; c_[59*max_+35] =    21631432489303455ull; c_[59*max_+36] =    14420954992868970ull; c_[59*max_+37] =     8964377427999630ull;
      c_[59*max_+38] =     5189902721473470ull; c_[59*max_+39] =     2794563003870330ull; c_[59*max_+40] =     1397281501935165ull; c_[59*max_+41] =      647520696018735ull;
      c_[59*max_+42] =      277508869722315ull; c_[59*max_+43] =      109712808959985ull; c_[59*max_+44] =       39895566894540ull; c_[59*max_+45] =       13298522298180ull;
      c_[59*max_+46] =        4047376351620ull; c_[59*max_+47] =        1119487075980ull; c_[59*max_+48] =         279871768995ull; c_[59*max_+49] =          62828356305ull;
      c_[59*max_+50] =          12565671261ull; c_[59*max_+51] =           2217471399ull; c_[59*max_+52] =            341149446ull; c_[59*max_+53] =             45057474ull;
      c_[59*max_+54] =              5006386ull; c_[59*max_+55] =               455126ull; c_[59*max_+56] =                32509ull; c_[59*max_+57] =                 1711ull;
      c_[59*max_+58] =                   59ull; c_[59*max_+59] =                    1ull; c_[60*max_+ 0] =                    1ull; c_[60*max_+ 1] =                   60ull;
      c_[60*max_+ 2] =                 1770ull; c_[60*max_+ 3] =                34220ull; c_[60*max_+ 4] =               487635ull; c_[60*max_+ 5] =              5461512ull;
      c_[60*max_+ 6] =             50063860ull; c_[60*max_+ 7] =            386206920ull; c_[60*max_+ 8] =           2558620845ull; c_[60*max_+ 9] =          14783142660ull;
      c_[60*max_+10] =          75394027566ull; c_[60*max_+11] =         342700125300ull; c_[60*max_+12] =        1399358844975ull; c_[60*max_+13] =        5166863427600ull;
      c_[60*max_+14] =       17345898649800ull; c_[60*max_+15] =       53194089192720ull; c_[60*max_+16] =      149608375854525ull; c_[60*max_+17] =      387221678682300ull;
      c_[60*max_+18] =      925029565741050ull; c_[60*max_+19] =     2044802197953900ull; c_[60*max_+20] =     4191844505805495ull; c_[60*max_+21] =     7984465725343800ull;
      c_[60*max_+22] =    14154280149473100ull; c_[60*max_+23] =    23385332420868600ull; c_[60*max_+24] =    36052387482172425ull; c_[60*max_+25] =    51915437974328292ull;
      c_[60*max_+26] =    69886166503903470ull; c_[60*max_+27] =    88004802264174740ull; c_[60*max_+28] =   103719945525634515ull; c_[60*max_+29] =   114449595062769120ull;
      c_[60*max_+30] =   118264581564861424ull; c_[60*max_+31] =   114449595062769120ull; c_[60*max_+32] =   103719945525634515ull; c_[60*max_+33] =    88004802264174740ull;
      c_[60*max_+34] =    69886166503903470ull; c_[60*max_+35] =    51915437974328292ull; c_[60*max_+36] =    36052387482172425ull; c_[60*max_+37] =    23385332420868600ull;
      c_[60*max_+38] =    14154280149473100ull; c_[60*max_+39] =     7984465725343800ull; c_[60*max_+40] =     4191844505805495ull; c_[60*max_+41] =     2044802197953900ull;
      c_[60*max_+42] =      925029565741050ull; c_[60*max_+43] =      387221678682300ull; c_[60*max_+44] =      149608375854525ull; c_[60*max_+45] =       53194089192720ull;
      c_[60*max_+46] =       17345898649800ull; c_[60*max_+47] =        5166863427600ull; c_[60*max_+48] =        1399358844975ull; c_[60*max_+49] =         342700125300ull;
      c_[60*max_+50] =          75394027566ull; c_[60*max_+51] =          14783142660ull; c_[60*max_+52] =           2558620845ull; c_[60*max_+53] =            386206920ull;
      c_[60*max_+54] =             50063860ull; c_[60*max_+55] =              5461512ull; c_[60*max_+56] =               487635ull; c_[60*max_+57] =                34220ull;
      c_[60*max_+58] =                 1770ull; c_[60*max_+59] =                   60ull; c_[60*max_+60] =                    1ull; c_[61*max_+ 0] =                    1ull;
      c_[61*max_+ 1] =                   61ull; c_[61*max_+ 2] =                 1830ull; c_[61*max_+ 3] =                35990ull; c_[61*max_+ 4] =               521855ull;
      c_[61*max_+ 5] =              5949147ull; c_[61*max_+ 6] =             55525372ull; c_[61*max_+ 7] =            436270780ull; c_[61*max_+ 8] =           2944827765ull;
      c_[61*max_+ 9] =          17341763505ull; c_[61*max_+10] =          90177170226ull; c_[61*max_+11] =         418094152866ull; c_[61*max_+12] =        1742058970275ull;
      c_[61*max_+13] =        6566222272575ull; c_[61*max_+14] =       22512762077400ull; c_[61*max_+15] =       70539987842520ull; c_[61*max_+16] =      202802465047245ull;
      c_[61*max_+17] =      536830054536825ull; c_[61*max_+18] =     1312251244423350ull; c_[61*max_+19] =     2969831763694950ull; c_[61*max_+20] =     6236646703759395ull;
      c_[61*max_+21] =    12176310231149295ull; c_[61*max_+22] =    22138745874816900ull; c_[61*max_+23] =    37539612570341700ull; c_[61*max_+24] =    59437719903041025ull;
      c_[61*max_+25] =    87967825456500717ull; c_[61*max_+26] =   121801604478231762ull; c_[61*max_+27] =   157890968768078210ull; c_[61*max_+28] =   191724747789809255ull;
      c_[61*max_+29] =   218169540588403635ull; c_[61*max_+30] =   232714176627630544ull; c_[61*max_+31] =   232714176627630544ull; c_[61*max_+32] =   218169540588403635ull;
      c_[61*max_+33] =   191724747789809255ull; c_[61*max_+34] =   157890968768078210ull; c_[61*max_+35] =   121801604478231762ull; c_[61*max_+36] =    87967825456500717ull;
      c_[61*max_+37] =    59437719903041025ull; c_[61*max_+38] =    37539612570341700ull; c_[61*max_+39] =    22138745874816900ull; c_[61*max_+40] =    12176310231149295ull;
      c_[61*max_+41] =     6236646703759395ull; c_[61*max_+42] =     2969831763694950ull; c_[61*max_+43] =     1312251244423350ull; c_[61*max_+44] =      536830054536825ull;
      c_[61*max_+45] =      202802465047245ull; c_[61*max_+46] =       70539987842520ull; c_[61*max_+47] =       22512762077400ull; c_[61*max_+48] =        6566222272575ull;
      c_[61*max_+49] =        1742058970275ull; c_[61*max_+50] =         418094152866ull; c_[61*max_+51] =          90177170226ull; c_[61*max_+52] =          17341763505ull;
      c_[61*max_+53] =           2944827765ull; c_[61*max_+54] =            436270780ull; c_[61*max_+55] =             55525372ull; c_[61*max_+56] =              5949147ull;
      c_[61*max_+57] =               521855ull; c_[61*max_+58] =                35990ull; c_[61*max_+59] =                 1830ull; c_[61*max_+60] =                   61ull;
      c_[61*max_+61] =                    1ull; c_[62*max_+ 0] =                    1ull; c_[62*max_+ 1] =                   62ull; c_[62*max_+ 2] =                 1891ull;
      c_[62*max_+ 3] =                37820ull; c_[62*max_+ 4] =               557845ull; c_[62*max_+ 5] =              6471002ull; c_[62*max_+ 6] =             61474519ull;
      c_[62*max_+ 7] =            491796152ull; c_[62*max_+ 8] =           3381098545ull; c_[62*max_+ 9] =          20286591270ull; c_[62*max_+10] =         107518933731ull;
      c_[62*max_+11] =         508271323092ull; c_[62*max_+12] =        2160153123141ull; c_[62*max_+13] =        8308281242850ull; c_[62*max_+14] =       29078984349975ull;
      c_[62*max_+15] =       93052749919920ull; c_[62*max_+16] =      273342452889765ull; c_[62*max_+17] =      739632519584070ull; c_[62*max_+18] =     1849081298960175ull;
      c_[62*max_+19] =     4282083008118300ull; c_[62*max_+20] =     9206478467454345ull; c_[62*max_+21] =    18412956934908690ull; c_[62*max_+22] =    34315056105966195ull;
      c_[62*max_+23] =    59678358445158600ull; c_[62*max_+24] =    96977332473382725ull; c_[62*max_+25] =   147405545359541742ull; c_[62*max_+26] =   209769429934732479ull;
      c_[62*max_+27] =   279692573246309972ull; c_[62*max_+28] =   349615716557887465ull; c_[62*max_+29] =   409894288378212890ull; c_[62*max_+30] =   450883717216034179ull;
      c_[62*max_+31] =   465428353255261088ull; c_[62*max_+32] =   450883717216034179ull; c_[62*max_+33] =   409894288378212890ull; c_[62*max_+34] =   349615716557887465ull;
      c_[62*max_+35] =   279692573246309972ull; c_[62*max_+36] =   209769429934732479ull; c_[62*max_+37] =   147405545359541742ull; c_[62*max_+38] =    96977332473382725ull;
      c_[62*max_+39] =    59678358445158600ull; c_[62*max_+40] =    34315056105966195ull; c_[62*max_+41] =    18412956934908690ull; c_[62*max_+42] =     9206478467454345ull;
      c_[62*max_+43] =     4282083008118300ull; c_[62*max_+44] =     1849081298960175ull; c_[62*max_+45] =      739632519584070ull; c_[62*max_+46] =      273342452889765ull;
      c_[62*max_+47] =       93052749919920ull; c_[62*max_+48] =       29078984349975ull; c_[62*max_+49] =        8308281242850ull; c_[62*max_+50] =        2160153123141ull;
      c_[62*max_+51] =         508271323092ull; c_[62*max_+52] =         107518933731ull; c_[62*max_+53] =          20286591270ull; c_[62*max_+54] =           3381098545ull;
      c_[62*max_+55] =            491796152ull; c_[62*max_+56] =             61474519ull; c_[62*max_+57] =              6471002ull; c_[62*max_+58] =               557845ull;
      c_[62*max_+59] =                37820ull; c_[62*max_+60] =                 1891ull; c_[62*max_+61] =                   62ull; c_[62*max_+62] =                    1ull;
      c_[63*max_+ 0] =                    1ull; c_[63*max_+ 1] =                   63ull; c_[63*max_+ 2] =                 1953ull; c_[63*max_+ 3] =                39711ull;
      c_[63*max_+ 4] =               595665ull; c_[63*max_+ 5] =              7028847ull; c_[63*max_+ 6] =             67945521ull; c_[63*max_+ 7] =            553270671ull;
      c_[63*max_+ 8] =           3872894697ull; c_[63*max_+ 9] =          23667689815ull; c_[63*max_+10] =         127805525001ull; c_[63*max_+11] =         615790256823ull;
      c_[63*max_+12] =        2668424446233ull; c_[63*max_+13] =       10468434365991ull; c_[63*max_+14] =       37387265592825ull; c_[63*max_+15] =      122131734269895ull;
      c_[63*max_+16] =      366395202809685ull; c_[63*max_+17] =     1012974972473835ull; c_[63*max_+18] =     2588713818544245ull; c_[63*max_+19] =     6131164307078475ull;
      c_[63*max_+20] =    13488561475572645ull; c_[63*max_+21] =    27619435402363035ull; c_[63*max_+22] =    52728013040874885ull; c_[63*max_+23] =    93993414551124795ull;
      c_[63*max_+24] =   156655690918541325ull; c_[63*max_+25] =   244382877832924467ull; c_[63*max_+26] =   357174975294274221ull; c_[63*max_+27] =   489462003181042451ull;
      c_[63*max_+28] =   629308289804197437ull; c_[63*max_+29] =   759510004936100355ull; c_[63*max_+30] =   860778005594247069ull; c_[63*max_+31] =   916312070471295267ull;
      c_[63*max_+32] =   916312070471295267ull; c_[63*max_+33] =   860778005594247069ull; c_[63*max_+34] =   759510004936100355ull; c_[63*max_+35] =   629308289804197437ull;
      c_[63*max_+36] =   489462003181042451ull; c_[63*max_+37] =   357174975294274221ull; c_[63*max_+38] =   244382877832924467ull; c_[63*max_+39] =   156655690918541325ull;
      c_[63*max_+40] =    93993414551124795ull; c_[63*max_+41] =    52728013040874885ull; c_[63*max_+42] =    27619435402363035ull; c_[63*max_+43] =    13488561475572645ull;
      c_[63*max_+44] =     6131164307078475ull; c_[63*max_+45] =     2588713818544245ull; c_[63*max_+46] =     1012974972473835ull; c_[63*max_+47] =      366395202809685ull;
      c_[63*max_+48] =      122131734269895ull; c_[63*max_+49] =       37387265592825ull; c_[63*max_+50] =       10468434365991ull; c_[63*max_+51] =        2668424446233ull;
      c_[63*max_+52] =         615790256823ull; c_[63*max_+53] =         127805525001ull; c_[63*max_+54] =          23667689815ull; c_[63*max_+55] =           3872894697ull;
      c_[63*max_+56] =            553270671ull; c_[63*max_+57] =             67945521ull; c_[63*max_+58] =              7028847ull; c_[63*max_+59] =               595665ull;
      c_[63*max_+60] =                39711ull; c_[63*max_+61] =                 1953ull; c_[63*max_+62] =                   63ull; c_[63*max_+63] =                    1ull;
      c_[64*max_+ 0] =                    1ull; c_[64*max_+ 1] =                   64ull; c_[64*max_+ 2] =                 2016ull; c_[64*max_+ 3] =                41664ull;
      c_[64*max_+ 4] =               635376ull; c_[64*max_+ 5] =              7624512ull; c_[64*max_+ 6] =             74974368ull; c_[64*max_+ 7] =            621216192ull;
      c_[64*max_+ 8] =           4426165368ull; c_[64*max_+ 9] =          27540584512ull; c_[64*max_+10] =         151473214816ull; c_[64*max_+11] =         743595781824ull;
      c_[64*max_+12] =        3284214703056ull; c_[64*max_+13] =       13136858812224ull; c_[64*max_+14] =       47855699958816ull; c_[64*max_+15] =      159518999862720ull;
      c_[64*max_+16] =      488526937079580ull; c_[64*max_+17] =     1379370175283520ull; c_[64*max_+18] =     3601688791018080ull; c_[64*max_+19] =     8719878125622720ull;
      c_[64*max_+20] =    19619725782651120ull; c_[64*max_+21] =    41107996877935680ull; c_[64*max_+22] =    80347448443237920ull; c_[64*max_+23] =   146721427591999680ull;
      c_[64*max_+24] =   250649105469666120ull; c_[64*max_+25] =   401038568751465792ull; c_[64*max_+26] =   601557853127198688ull; c_[64*max_+27] =   846636978475316672ull;
      c_[64*max_+28] =  1118770292985239888ull; c_[64*max_+29] =  1388818294740297792ull; c_[64*max_+30] =  1620288010530347424ull; c_[64*max_+31] =  1777090076065542336ull;
      c_[64*max_+32] =  1832624140942590534ull; c_[64*max_+33] =  1777090076065542336ull; c_[64*max_+34] =  1620288010530347424ull; c_[64*max_+35] =  1388818294740297792ull;
      c_[64*max_+36] =  1118770292985239888ull; c_[64*max_+37] =   846636978475316672ull; c_[64*max_+38] =   601557853127198688ull; c_[64*max_+39] =   401038568751465792ull;
      c_[64*max_+40] =   250649105469666120ull; c_[64*max_+41] =   146721427591999680ull; c_[64*max_+42] =    80347448443237920ull; c_[64*max_+43] =    41107996877935680ull;
      c_[64*max_+44] =    19619725782651120ull; c_[64*max_+45] =     8719878125622720ull; c_[64*max_+46] =     3601688791018080ull; c_[64*max_+47] =     1379370175283520ull;
      c_[64*max_+48] =      488526937079580ull; c_[64*max_+49] =      159518999862720ull; c_[64*max_+50] =       47855699958816ull; c_[64*max_+51] =       13136858812224ull;
      c_[64*max_+52] =        3284214703056ull; c_[64*max_+53] =         743595781824ull; c_[64*max_+54] =         151473214816ull; c_[64*max_+55] =          27540584512ull;
      c_[64*max_+56] =           4426165368ull; c_[64*max_+57] =            621216192ull; c_[64*max_+58] =             74974368ull; c_[64*max_+59] =              7624512ull;
      c_[64*max_+60] =               635376ull; c_[64*max_+61] =                41664ull; c_[64*max_+62] =                 2016ull; c_[64*max_+63] =                   64ull;
      c_[64*max_+64] =                    1ull;
    }
    size_t operator()(const int i, const int j) const { assert(i >= 0 && j >= 0); return c_[i*max_+j]; }
};

}

#endif

#if 0
// Code to generate this table
#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpreal.h"
using namespace mpfr;
using namespace std;

mpreal fac(const mpreal i) {
  if (i < 1.1) { return 1; }
  else { return i * fac(i-1); }
}

int main() {
  mpreal::set_default_prec(10000);
  for (int i=0, ij=0; i!=65; ++i) {
    for (int j=0; j!=i+1; ++j, ++ij) {
      mpreal ii(i);
      mpreal jj(j);
      mpreal kk = ii - jj;
      mpreal out = fac(ii)/fac(jj)/fac(kk);
      // out should be an integer
      if (abs(out-mpreal(out.toULong())) > 1.0e-10) assert(false);
      cout << "c_[" << setw(2) << i << "*max_+" << setw(2) << j << "] = " << setw(20) << out.toULong() << "ull;";
      if (ij % 4 == 3) cout << endl;
      else cout << " ";
    }
  }
  return 0;
}
#endif
