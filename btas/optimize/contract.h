// TODO Only ColMajor so far.
#ifndef BTAS_OPTIMIZE_CONTRACT_H
#define BTAS_OPTIMIZE_CONTRACT_H

#include <sstream>
#include <stdexcept>
#include <cassert>
#include <btas/generic/gemm_impl.h>
#include <btas/generic/scal_impl.h>

namespace btas {

template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_211(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC, const bool conjgA, const bool conjgB) {
  assert(aA.size() == 2 && aB.size() == 1 && aC.size() == 1);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());
  if (conjgB) throw std::logic_error("complex conjugation of 1-index tensors is not considered in contract_211");

  const bool notrans = aB[0] == aA[1];
  if (notrans && conjgA) throw std::logic_error("contract_211 not sure what to do");

  auto cA = notrans ? CblasNoTrans : (conjgA ? CblasConjTrans : CblasTrans);
  assert((notrans && aA[0] == aC[0]) || (aB[0] == aA[0] && aA[1] == aC[0]));
  gemv_impl<true>::call(CblasColMajor, cA, A.extent(0), A.extent(1), 1.0, &*A.begin(), A.extent(0), &*B.begin(), 1, 0.0, &*C.begin(), 1);
}


template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_222(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC, const bool conjgA, const bool conjgB) {
  // TODO we do not consider complex matrixces yet.
  assert(aA.size() == 2 && aB.size() == 2 && aC.size() == 2);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());
  if (std::find(aA.begin(), aA.end(), aC.front()) != aA.end()) {
    // then multiply A * B -> C
    const bool notransA = aA.front() == aC.front();
    const bool notransB = aB.front() != aC.back();
    if (notransA && conjgA) throw std::logic_error("contract_211 not sure what to do (A)");
    if (notransB && conjgB) throw std::logic_error("contract_211 not sure what to do (B)");

    const auto   cA     = notransA ? CblasNoTrans : (conjgA ? CblasConjTrans : CblasTrans);
    const size_t condim = notransA ? A.extent(1) : A.extent(0);
    assert(std::find(aB.begin(), aB.end(), aC.back()) != aB.end());

    const auto   cB     = notransB ? CblasNoTrans : (conjgB ? CblasConjTrans : CblasTrans);
    assert((notransA ? aA.back() : aA.front()) == (notransB ? aB.front() : aB.back()));

    gemm_impl<true>::call(CblasColMajor, cA, cB, C.extent(0), C.extent(1), condim,
                          alpha, &*A.begin(), A.extent(0), &*B.begin(), B.extent(0), beta, &*C.begin(), C.extent(0));
  } else {
    contract_222(alpha, B, aB, A, aA, beta, C, aC, conjgB, conjgA);
  }
}


template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_323(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC, const bool conjgA, const bool conjgB) {
  assert(aA.size() == 3 && aB.size() == 2 && aC.size() == 3);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());
  if (conjgA) throw std::logic_error("complex conjugation of 3-index tensors is not considered in contract_323");

  // TODO this function is limited to special cases where one of three indices of A will be replaced in C. Permuation is not considered so far.
  // first idenfity which indices to be rotated
  int irot = -1;
  for (int i = 0; i != 3; ++i)
    if (aA[i] != aC[i]) {
      assert(irot < 0);
      irot = i;
    } else
      assert(A.extent(i) == C.extent(i));

  if (irot == 0) {
    // in this case multiply from front
    const bool notrans = aB.back() == aA.front();
    assert(notrans || aB.front() == aA.front());
    if (notrans && conjgB) throw std::logic_error("contract_323 irot = 0 not sure what to do");

    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : (conjgB ? CblasConjTrans : CblasTrans);
    assert(notrans ? B.extent(1) : B.extent(0) == A.extent(0));

    gemm_impl<true>::call(CblasColMajor, cB, cA, C.extent(0), C.extent(1)*C.extent(2), A.extent(0),
                          alpha, &*B.begin(), B.extent(0), &*A.begin(), A.extent(0), beta, &*C.begin(), C.extent(0));
  } else if (irot == 1) {
    // in this case we loop over the last index of A
    const bool notrans = aB.front() == aA[1];
    assert(notrans || aB.back() == aA[1]);
    if (notrans && conjgB) throw std::logic_error("contract_323 irot = 1 not sure what to do");

    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : (conjgB ? CblasConjTrans : CblasTrans);
    assert(notrans ? B.extent(0) : B.extent(1) == A.extent(1));
    const size_t ablock = A.extent(0)*A.extent(1);
    const size_t cblock = C.extent(0)*C.extent(1);

    for (int i = 0; i != A.extent(2); ++i)
      gemm_impl<true>::call(CblasColMajor, cA, cB, C.extent(0), C.extent(1), A.extent(1),
                            alpha, &*A.begin()+i*ablock, A.extent(0), &*B.begin(), B.extent(0), beta, &*C.begin()+i*cblock, C.extent(0));
  } else if (irot == 2) {
    // in this case multiply from back
    const bool notrans = aB.front() == aA[2];
    assert(notrans || aB.back() == aA[2]);
    if (notrans && conjgB) throw std::logic_error("contract_323 irot = 2 not sure what to do");

    const auto cA = CblasNoTrans;
    const auto cB = notrans ? CblasNoTrans : (conjgB ? CblasConjTrans : CblasTrans);
    assert(notrans ? B.extent(0) : B.extent(1) == A.extent(2));
    gemm_impl<true>::call(CblasColMajor, cA, cB, C.extent(0)*C.extent(1), C.extent(2), A.extent(2),
                          alpha, &*A.begin(), A.extent(0)*A.extent(1), &*B.begin(), B.extent(0), beta, &*C.begin(), C.extent(0)*C.extent(1));
  } else {
    assert(false);
  }
}


template<typename _T, class _TensorA, class _TensorB, class _TensorC,
         typename _UA, typename _UB, typename _UC
        >
void contract_332(const _T& alpha, const _TensorA& A, const btas::varray<_UA>& aA, const _TensorB& B, const btas::varray<_UB>& aB,
                  const _T& beta, _TensorC& C, const btas::varray<_UC>& aC, const bool conjgA, const bool conjgB) {
  assert(aA.size() == 3 && aB.size() == 3 && aC.size() == 2);
  assert(A.range().ordinal().contiguous() && B.range().ordinal().contiguous() && C.range().ordinal().contiguous());

  const bool back2  = aA[0] == aB[0] && aA[1] == aB[1];
  const bool front2 = aA[1] == aB[1] && aA[2] == aB[2];
  const bool mid2   = aA[0] == aB[0] && aA[2] == aB[2];
  if (back2) {
    const bool swap = aC[0] == aB[2];
    assert(swap || aC[0] == aA[2]);
    if (!swap) {
      assert(A.extent(0)*A.extent(1) == B.extent(0)*B.extent(1) && A.extent(2) == C.extent(0) && B.extent(2) == C.extent(1));
      assert(!conjgB);
      gemm_impl<true>::call(CblasColMajor, conjgA ? CblasConjTrans : CblasTrans, CblasNoTrans, C.extent(0), C.extent(1), A.extent(0)*A.extent(1),
                            alpha, &*A.begin(), A.extent(0)*A.extent(1), &*B.begin(), B.extent(0)*B.extent(1), beta, &*C.begin(), C.extent(0));
    } else {
      assert(A.extent(0)*A.extent(1) == B.extent(0)*B.extent(1) && B.extent(2) == C.extent(0) && A.extent(2) == C.extent(1));
      assert(!conjgA);
      gemm_impl<true>::call(CblasColMajor, conjgB ? CblasConjTrans : CblasTrans, CblasNoTrans, C.extent(0), C.extent(1), A.extent(0)*A.extent(1),
                            alpha, &*B.begin(), B.extent(0)*B.extent(1), &*A.begin(), A.extent(0)*A.extent(1), beta, &*C.begin(), C.extent(0));
    }
  } else if (front2) {
    const bool swap = aC[0] == aB[0];
    assert(swap || aC[0] == aA[0]);
    if (!swap) {
      assert(A.extent(1)*A.extent(2) == B.extent(1)*B.extent(2) && A.extent(0) == C.extent(0) && B.extent(0) == C.extent(1));
      assert(!conjgA);
      gemm_impl<true>::call(CblasColMajor, CblasNoTrans, conjgB ? CblasConjTrans : CblasTrans, C.extent(0), C.extent(1), A.extent(1)*A.extent(2),
                            alpha, &*A.begin(), A.extent(0), &*B.begin(), B.extent(0), beta, &*C.begin(), C.extent(0));
    } else {
      assert(A.extent(1)*A.extent(2) == B.extent(1)*B.extent(2) && B.extent(0) == C.extent(0) && A.extent(0) == C.extent(1));
      assert(!conjgB);
      gemm_impl<true>::call(CblasColMajor, CblasNoTrans, conjgA ? CblasConjTrans : CblasTrans, C.extent(0), C.extent(1), A.extent(1)*A.extent(2),
                            alpha, &*B.begin(), B.extent(0), &*A.begin(), A.extent(0), beta, &*C.begin(), C.extent(0));
    }
  } else if (mid2) {
    const bool swap = aC[0] == aB[1];
    assert(swap || aC[0] == aA[1]);
    const size_t ablock = A.extent(0)*A.extent(1);
    const size_t bblock = B.extent(0)*B.extent(1);
    scal(C.size(), beta, &*C.begin(), 1);
    if (!swap) {
      assert(A.extent(0) == B.extent(0) && A.extent(2) == B.extent(2) && A.extent(1) == C.extent(0) && B.extent(1) == C.extent(1));
      assert(!conjgB);
      for (int i = 0; i != A.extent(2); ++i)
        gemm_impl<true>::call(CblasColMajor, conjgA ? CblasConjTrans : CblasTrans, CblasNoTrans, C.extent(0), C.extent(1), A.extent(0),
                              alpha, &*A.begin()+i*ablock, A.extent(0), &*B.begin()+i*bblock, B.extent(0), 1.0, &*C.begin(), C.extent(0));
    } else {
      assert(A.extent(0) == B.extent(0) && A.extent(2) == B.extent(2) && B.extent(1) == C.extent(0) && A.extent(1) == C.extent(1));
      assert(!conjgA);
      for (int i = 0; i != A.extent(2); ++i)
        gemm_impl<true>::call(CblasColMajor, conjgB ? CblasConjTrans : CblasTrans, CblasNoTrans, C.extent(0), C.extent(1), A.extent(0),
                              alpha, &*B.begin()+i*bblock, B.extent(0), &*A.begin()+i*ablock, A.extent(0), 1.0, &*C.begin(), C.extent(0));
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
    (_TensorA::range_type::order == CblasColMajor) & //checking if A, B, and C are all Colomn major
    (_TensorB::range_type::order == CblasColMajor) & //checking if A, B, and C are all Colomn major
    (_TensorC::range_type::order == CblasColMajor) & //checking if A, B, and C are all Colomn major
    std::is_same<typename std::remove_cv<typename _TensorA::value_type>::type, typename std::remove_cv<typename _TensorB::value_type>::type>::value &
    std::is_same<typename std::remove_cv<typename _TensorA::value_type>::type, typename std::remove_cv<typename _TensorC::value_type>::type>::value &
    (std::is_same<typename std::remove_cv<typename _TensorA::value_type>::type, double>::value
  or std::is_same<typename std::remove_cv<typename _TensorA::value_type>::type, std::complex<double>>::value)
  >::type
>
void contract(
  const _T& alpha,
  const _TensorA& A, std::initializer_list<_UA> aA,
  const _TensorB& B, std::initializer_list<_UB> aB,
  const _T& beta,
        _TensorC& C, std::initializer_list<_UC> aC,
  const bool conjgA = false,
  const bool conjgB = false) {

  assert(A.rank() == aA.size());
  assert(B.rank() == aB.size());
  assert(C.rank() == aC.size());

  if (A.rank() == 2 && B.rank() == 1 && C.rank() == 1) {
    contract_211(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC), conjgA, conjgB);
  } else if (A.rank() == 1 && B.rank() == 2 && C.rank() == 1) {
    contract_211(alpha, B, btas::varray<_UA>(aB), A, btas::varray<_UB>(aA), beta, C, btas::varray<_UC>(aC), conjgB, conjgA);
  } else if (A.rank() == 2 && B.rank() == 2 && C.rank() == 2) {
    contract_222(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC), conjgA, conjgB);
  } else if (A.rank() == 3 && B.rank() == 2 && C.rank() == 3) {
    contract_323(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC), conjgA, conjgB);
  } else if (A.rank() == 2 && B.rank() == 3 && C.rank() == 3) {
    contract_323(alpha, B, btas::varray<_UA>(aB), A, btas::varray<_UB>(aA), beta, C, btas::varray<_UC>(aC), conjgB, conjgA);
  } else if (A.rank() == 3 && B.rank() == 3 && C.rank() == 2) {
    contract_332(alpha, A, btas::varray<_UA>(aA), B, btas::varray<_UB>(aB), beta, C, btas::varray<_UC>(aC), conjgA, conjgB);
  } else {
    std::stringstream ss; ss << "not yet implemented: rank(A): " << A.rank() << " rank(B): " << B.rank() << " rank(C): " << C.rank();
    throw std::logic_error(ss.str());
  }
}

} //namespace btas

#endif
