//
// Author: Toru Shiozaki
// Date  : April 2009
// 

#include "angular_index.h"
#include <cassert>

using namespace std;
using namespace boost;

const string Angular_Index::show() const {
  string out;

  for (int i = 0; i != get<0>(index_); ++i) out += "x";
  for (int i = 0; i != get<1>(index_); ++i) out += "y";
  for (int i = 0; i != get<2>(index_); ++i) out += "z";
  if (out.empty()) out += "0";

  return out;
}


tuple<Angular_Pair, Angular_Pair, int> Angular_Pair::hrr_formula() const {

  Angular_Index mine1 = indices_.first;
  Angular_Index mine2 = indices_.second;

  if (mine2.x() != 0) {
    // for first term
    Angular_Index term11(mine1.x() + 1, mine1.y(), mine1.z());  
    Angular_Index term12(mine2.x() - 1, mine2.y(), mine2.z());  
    Angular_Pair term1(make_pair(term11, term12));
    // for second term
    Angular_Index term21(mine1.x(), mine1.y(), mine1.z());  
    Angular_Index term22(mine2.x() - 1, mine2.y(), mine2.z());  
    Angular_Pair term2(make_pair(term21, term22));

    return make_tuple(term1, term2, 0);
  } else if (mine2.y() != 0) {
    // for first term
    Angular_Index term11(mine1.x(), mine1.y() + 1, mine1.z());  
    Angular_Index term12(mine2.x(), mine2.y() - 1, mine2.z());  
    Angular_Pair term1(make_pair(term11, term12));
    // for second term
    Angular_Index term21(mine1.x(), mine1.y(), mine1.z());  
    Angular_Index term22(mine2.x(), mine2.y() - 1, mine2.z());  
    Angular_Pair term2(make_pair(term21, term22));

    return make_tuple(term1, term2, 1);
  } else {
    assert(mine2.z() != 0);
    // for first term
    Angular_Index term11(mine1.x(), mine1.y(), mine1.z() + 1);  
    Angular_Index term12(mine2.x(), mine2.y(), mine2.z() - 1);  
    Angular_Pair term1(make_pair(term11, term12));
    // for second term
    Angular_Index term21(mine1.x(), mine1.y(), mine1.z());  
    Angular_Index term22(mine2.x(), mine2.y(), mine2.z() - 1);  
    Angular_Pair term2(make_pair(term21, term22));

    return make_tuple(term1, term2, 2);
  }

}


const string Angular_Pair::show() const {
  string out;
  out += "a" + indices_.first.show() + "_" + indices_.second.show();
  return out;
}
