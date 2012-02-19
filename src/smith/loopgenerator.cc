//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include <src/smith/loopgenerator.h>
#include <algorithm>

using namespace SMITH;
using namespace std;


vector<vector<Index> > LoopGenerator::block_loop() const {
  // first, make a status vector 
  std::vector<int> stat(loop_.size());
  std::vector<int> max(loop_.size());
  {
    auto j = loop_.begin();
    for (auto i = max.begin(); i != max.end(); ++i, ++j) *i = j->nblock(); 
  }

  vector<vector<Index> > out;

  do {
    vector<Index> tmp;
    auto l = loop_.begin();
    for (auto k = stat.begin(); k != stat.end(); ++k, ++l)
      tmp.push_back(l->range(*k)); 
    out.push_back(tmp);

    auto j = stat.begin();
    auto i = max.begin();
    while (j != stat.end() && (++*j) == *i) {
      *j = 0; 
      ++j;
      ++i;
    }
  } while (*max_element(stat.begin(), stat.end()) > 0);

  return out;
}
