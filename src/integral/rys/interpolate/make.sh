//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
#g++ -O0 -c ../_breitroot_*.cc
#g++ -O0 -c ../_spin2root_*.cc

#g++ -std=c++11 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o -o gen
#g++ -Ofast -std=c++11 -DBREIT -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen
#g++ -Ofast -std=c++11 -DSPIN2 -lblas -llapack -lgmp -lmpfr *.cc *.c _breit*.o _spin2*.o -o gen

# for debug
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