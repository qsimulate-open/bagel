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
  // first compute "N" and "T" things
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
  } else {
    throw std::logic_error("not yet implemented");
  }
}

} //namespace btas

#endif
