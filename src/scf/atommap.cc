//
// Author: Toru Shiozaki
// May 2009
//

#include <src/scf/atommap.h>

using namespace std;

AtomMap::AtomMap () {
  atommap.insert(make_pair("H" ,  1));
  atommap.insert(make_pair("He",  2));
  atommap.insert(make_pair("Li",  3));
  atommap.insert(make_pair("Be",  4));
  atommap.insert(make_pair("B" ,  5));
  atommap.insert(make_pair("C" ,  6));
  atommap.insert(make_pair("N" ,  7));
  atommap.insert(make_pair("O" ,  8));
  atommap.insert(make_pair("F" ,  9));
  atommap.insert(make_pair("Ne", 10));

  angmap.insert(make_pair("s", 0));
  angmap.insert(make_pair("p", 1));
  angmap.insert(make_pair("d", 2));
  angmap.insert(make_pair("f", 3));
  angmap.insert(make_pair("g", 4));
  angmap.insert(make_pair("h", 5));
  angmap.insert(make_pair("i", 6));
  angmap.insert(make_pair("j", 7));
  angmap.insert(make_pair("k", 8));
  angmap.insert(make_pair("l", 9));

}


AtomMap::~AtomMap () {

}
