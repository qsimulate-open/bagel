//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <string>
#include <cmath>

class VRR {
  private:
    // target angular momentum
    int a_, c_;

    // rank of Rys quadruture
    int rank_;

    // generating functions for each case
    const std::pair<std::string, int> vrr00 (         ) const; 
    const std::pair<std::string, int> vrrn0 (const int) const; 
    const std::pair<std::string, int> vrr0m (const int) const; 
    const std::pair<std::string, int> vrr11 (         ) const; 
    const std::pair<std::string, int> vrrn1 (const int) const;
    const std::pair<std::string, int> vrr1m (const int) const;
    const std::pair<std::string, int> vrrnm (const int, const int) const;

  public:
    // CAUTION : a change due to gradient (since _i and _j are not used at the same time). 
    VRR(const int _i, const int _j): a_(_i), c_(_j), rank_(::ceil(0.5 * (_i + _j))) { }; 
    //VRR(const int _i, const int _j): a_(_i), c_(_j), rank_(::ceil(0.5 * (_i + _j + 1))) { }; 
    ~VRR() { };

    const std::string dump(const std::string) const;

};

