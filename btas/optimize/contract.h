// I implement these things for BAGEL first, and will make it generic
#ifndef BTAS_OPTIMIZE_CONTRACT_H
#define BTAS_OPTIMIZE_CONTRACT_H

#include <stdexcept>
#include <cassert>
#include <btas/generic/gemm_impl.h>

namespace btas {

template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_222(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC) {
  // TODO we do not consider complex matrixces yet.
  assert(aA.size() == 2 && aB.size() == 2 && aC.size() == 2);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());
  if (std::find(aA.begin(), aA.end(), aC.front()) != aA.end()) {
    // then multiply A * B -> C
    const bool notrans  = aA.front() == aC.front();
    const auto   cA     = notrans ? CblasNoTrans : CblasTrans;
    const size_t condim = notrans ? A.range(1).size() : A.range(0).size();
    assert(std::find(aB.begin(), aB.end(), aC.back()) != aB.end());
    const auto   cB     = aB.front() == aC.back() ? CblasTrans : CblasNoTrans;
    assert((cA == CblasNoTrans ? aA.back() : aA.front()) == (cB == CblasTrans ? aB.back() : aB.front()));

    gemm_impl<true>::call(CblasColMajor, cA, cB, C.range(0).size(), C.range(1).size(), condim,
                          alpha, A.data(), A.range(0).size(), B.data(), B.range(0).size(), beta, C.data(), C.range(0).size());
  } else {
    contract_222(alpha, B, aB, A, aA, beta, C, aC);
  }
}


template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_323(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC) {
  assert(aA.size() == 3 && aB.size() == 2 && aC.size() == 3);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());

  // TODO this function is limited to special cases where one of three indices of A will be replaced in C. Permuation is not considered so far.
  // first idenfity which indices to be rotated
  int irot = -1;
  for (int i = 0; i != 3; ++i)
    if (aA[i] != aC[i]) {
      assert(irot < 0);
      irot = i;
    } else
      assert(A.range(i).size() == C.range(i).size());

  if (irot == 0) {
    // in this case multiply from front
    const bool notrans = aB.back() == aA.front();
    assert(notrans || aB.front() == aA.front());
    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : CblasTrans;
    assert(notrans ? B.range(1).size() : B.range(0).size() == A.range(0).size());

    gemm_impl<true>::call(CblasColMajor, cB, cA, C.range(0).size(), C.range(1).size()*C.range(2).size(), A.range(0).size(),
                          alpha, B.data(), B.range(0).size(), A.data(), A.range(0).size(), beta, C.data(), C.range(0).size());
  } else if (irot == 1) {
    // in this case we loop over the last index of A
    const bool notrans = aB.front() == aA[1];
    assert(notrans || aB.back() == aA[1]);
    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : CblasTrans;
    assert(notrans ? B.range(0).size() : B.range(1).size() == A.range(1).size());
    const size_t ablock = A.range(0).size()*A.range(1).size();
    const size_t cblock = C.range(0).size()*C.range(1).size();

    for (int i = 0; i != A.range(2).size(); ++i)
      gemm_impl<true>::call(CblasColMajor, cA, cB, C.range(0).size(), C.range(1).size(), A.range(1).size(),
                            alpha, A.data()+i*ablock, A.range(0).size(), B.data(), B.range(0).size(), beta, C.data()+i*cblock, C.range(0).size());
  } else if (irot == 2) {
    // in this case multiply from back
    const bool notrans = aB.front() == aA[2];
    assert(notrans || aB.back() == aA[2]);
    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : CblasTrans;
    assert(notrans ? B.range(0).size() : B.range(1).size() == A.range(2).size());
    gemm_impl<true>::call(CblasColMajor, cA, cB, C.range(0).size()*C.range(1).size(), C.range(2).size(), A.range(2).size(),
                          alpha, A.data(), A.range(0).size()*A.range(1).size(), B.data(), B.range(0).size(), beta, C.data(), C.range(0).size()*C.range(1).size());
  } else {
    assert(false);
  }
}


template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_332(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC) {
  assert(aA.size() == 3 && aB.size() == 3 && aC.size() == 2);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());

  const bool back2  = aA[0] == aB[0] && aA[1] == aB[1];
  const bool front2 = aA[1] == aB[1] && aA[2] == aB[2];
  const bool mid2   = aA[0] == aB[0] && aA[2] == aB[2];
  if (back2) {
    const bool swap = aC[0] == aB[2];
    assert(swap || aC[0] == aA[2]);
    if (!swap) {
      assert(A.range(0).size()*A.range(1).size() == B.range(0).size()*B.range(1).size() && A.range(2).size() == C.range(0).size() && B.range(2).size() == C.range(1).size());
      gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, C.range(0).size(), C.range(1).size(), A.range(0).size()*A.range(1).size(),
                            alpha, A.data(), A.range(0).size()*A.range(1).size(), B.data(), B.range(0).size()*B.range(1).size(), beta, C.data(), C.range(0).size());
    } else {
      assert(A.range(0).size()*A.range(1).size() == B.range(0).size()*B.range(1).size() && B.range(2).size() == C.range(0).size() && A.range(2).size() == C.range(1).size());
      gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, C.range(0).size(), C.range(1).size(), A.range(0).size()*A.range(1).size(),
                            alpha, B.data(), B.range(0).size()*B.range(1).size(), A.data(), A.range(0).size()*A.range(1).size(), beta, C.data(), C.range(0).size());
    }
  } else if (front2) {
    const bool swap = aC[0] == aB[0];
    assert(swap || aC[0] == aA[0]);
    if (!swap) {
      assert(A.range(1).size()*A.range(2).size() == B.range(1).size()*B.range(2).size() && A.range(0).size() == C.range(0).size() && B.range(0).size() == C.range(1).size());
      gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, C.range(0).size(), C.range(1).size(), A.range(1).size()*A.range(2).size(),
                            alpha, A.data(), A.range(0).size(), B.data(), B.range(0).size(), beta, C.data(), C.range(0).size());
    } else {
      assert(A.range(1).size()*A.range(2).size() == B.range(1).size()*B.range(2).size() && B.range(0).size() == C.range(0).size() && A.range(0).size() == C.range(1).size());
      gemm_impl<true>::call(CblasColMajor, CblasNoTrans, CblasTrans, C.range(0).size(), C.range(1).size(), A.range(1).size()*A.range(2).size(),
                            alpha, B.data(), B.range(0).size(), A.data(), A.range(0).size(), beta, C.data(), C.range(0).size());
    }
  } else if (mid2) {
    const bool swap = aC[0] == aB[1];
    assert(swap || aC[0] == aA[1]);
    const size_t ablock = A.range(0).size()*A.range(1).size();
    const size_t bblock = B.range(0).size()*B.range(1).size();
    C.scale(beta);
    if (!swap) {
      assert(A.range(0).size() == B.range(0).size() && A.range(2).size() == B.range(2).size() && A.range(1).size() == C.range(0).size() && B.range(1).size() == C.range(1).size());
      for (int i = 0; i != A.range(2).size(); ++i)
        gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, C.range(0).size(), C.range(1).size(), A.range(0).size(),
                              alpha, A.data()+i*ablock, A.range(0).size(), B.data()+i*bblock, B.range(0).size(), 1.0, C.data(), C.range(0).size());
    } else {
      assert(A.range(0).size() == B.range(0).size() && A.range(2).size() == B.range(2).size() && B.range(1).size() == C.range(0).size() && A.range(1).size() == C.range(1).size());
      for (int i = 0; i != A.range(2).size(); ++i)
        gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, C.range(0).size(), C.range(1).size(), A.range(0).size(),
                              alpha, B.data()+i*bblock, B.range(0).size(), A.data()+i*ablock, A.range(0).size(), 1.0, C.data(), C.range(0).size());
    }
  } else
    throw std::logic_error("not yet implemented");
}


template<
  typename _T,
  class _TensorA, class _TensorB, class _TensorC,
  typename _UA, typename _UB, typename _UC,
  class = typename std::enable_if<
    is_tensor<_TensorA>::value &
    is_tensor<_TensorB>::value &
    is_tensor<_TensorC>::value &
    _TensorA::range_type::order == CblasColMajor & //checking if A, B, and C are all Colomn major
    _TensorB::range_type::order == CblasColMajor & //checking if A, B, and C are all Colomn major
    _TensorC::range_type::order == CblasColMajor & //checking if A, B, and C are all Colomn major
    std::is_same<typename _TensorA::value_type, typename _TensorB::value_type>::value &
    std::is_same<typename _TensorA::value_type, typename _TensorC::value_type>::value
  >::type
>
void contract(
  const _T& alpha,
  const _TensorA& A, std::initializer_list<_UA> aA,
  const _TensorB& B, std::initializer_list<_UB> aB,
  const _T& beta,
        _TensorC& C, std::initializer_list<_UC> aC) {

  if (A.rank() == 2 && B.rank() == 2 && C.rank() == 2) {
    contract_222(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC));
  } else if (A.rank() == 3 && B.rank() == 2 && C.rank() == 3) {
    contract_323(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC));
  } else if (A.rank() == 2 && B.rank() == 3 && C.rank() == 3) {
    contract_323(alpha, B, btas::varray<_UA>(aB), A, btas::varray<_UB>(aA), beta, C, btas::varray<_UC>(aC));
  } else if (A.rank() == 3 && B.rank() == 3 && C.rank() == 2) {
    contract_332(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC));
  } else {
    throw std::logic_error("not yet implemented");
  }
}

} //namespace btas

#endif
