/*
This header will define templates allowing us to add and multiply polynomials.
It will be used to allow my ERI-testing file to compute data for nonzero angular momentum.
*/

#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

namespace ryan {

template <typename T>
class polynomial {

  protected:

  public:

    int rank;
    std::vector<T> coeff;

    // constructor
    polynomial (std::vector<T> coefficients) {
      rank = coefficients.size() - 1;
      coeff = coefficients;
    }

    // to print out polynomial
    void show () const {
      for (int i=rank; i>=0; i--) {
        std::cout << coeff[i];
        if (i==0) std::cout << std::endl;
        else if (i==1) std::cout << "x" << " + ";
        else std::cout << "x^" << i << " + ";
      }
    }

    // to plug in a particular value of the independent variable (x)
    T evaluate (T x) const {
      T out = 0.0;
      for (int i=rank; i>=0; i--) {
        out += coeff[i]*std::pow(x,i);
      }
      return out;
    }

};

template <typename T>
polynomial<T> multiply_polynomials (polynomial<T> polyA, polynomial<T> polyB) {

  assert (polyA.rank + 1 == polyA.coeff.size());
  assert (polyB.rank + 1 == polyB.coeff.size());
  const int rankP = polyA.rank + polyB.rank;

  std::vector<T> new_coeffs (rankP+1, 0.0);
  for (int j=0; j<=polyA.rank; j++) {
    for (int k=0; k<=polyB.rank; k++) {
      new_coeffs[j+k] += (polyA.coeff[j]*polyB.coeff[k]);
    }
  }

  polynomial<T> out (new_coeffs);
  return out;
}

template <typename T>
polynomial<T> scalar_polynomial (polynomial<T> polyA, T scalar) {
  std::vector<T> new_coeffs (polyA.rank + 1, 0.0);
  for (int j=0; j<=polyA.rank; j++) {
    new_coeffs[j] = scalar*polyA.coeff[j];
  }

  // Reduce rank if coefficient is zero
  const T zero = 0.0;
  while (new_coeffs.back() == zero && new_coeffs.size() > 1) new_coeffs.pop_back();

  polynomial<T> out (new_coeffs);
  return out;
}

template <typename T>
polynomial<T> add_polynomials (polynomial<T> polyA, polynomial<T> polyB) {

  const int rankS = std::max(polyA.rank,polyB.rank);
  const int lowrank = std::min(polyA.rank,polyB.rank);
  const T zero = 0.0;

  std::vector<T> new_coeffs (rankS+1, zero);
  for (int j=0; j<=lowrank; j++) {
    new_coeffs[j] = (polyA.coeff[j]+polyB.coeff[j]);
  }

  if (rankS > polyB.rank) {
    for (int j=lowrank+1; j<=rankS; j++) {
      new_coeffs[j] = polyA.coeff[j];
    }
  }

  if (rankS > polyA.rank) {
    for (int j=lowrank+1; j<=rankS; j++) {
      new_coeffs[j] = polyB.coeff[j];
    }
  }

  // Reduce rank if coefficient is zero
  while (new_coeffs.back() == zero && new_coeffs.size() > 1) new_coeffs.pop_back();

  polynomial<T> out (new_coeffs);
  return out;
}

template <typename T>
polynomial<T> subtract_polynomials (polynomial<T> polyA, polynomial<T> polyB) {
  T minus = -1.0;
  polynomial<T> minusB = scalar_polynomial (polyB, minus);
  return add_polynomials (polyA, minusB);
}


}

#endif
