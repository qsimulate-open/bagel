#ifndef __BTAS_CONTRACT_H
#define __BTAS_CONTRACT_H

#include <algorithm>
#include <memory>

#include <btas/types.h>
#include <btas/tensor_traits.h>

#include <btas/util/resize.h>

#include <btas/generic/scal_impl.h>
#include <btas/generic/gemv_impl.h>
#include <btas/generic/ger_impl.h>
#include <btas/generic/gemm_impl.h>
#include <btas/generic/permute.h>

namespace btas {

/// contract tensors; for example, Cijk = \sum_{m,p} Aimp * Bmjpk
///
/// Synopsis:
/// enum {j,k,l,m,n,o};
///
/// contract(alpha,A,{m,o,k,n},B,{l,k,j},C,beta,{l,n,m,o,j});
///
///       o       j           o j
///       |       |           | |
///   m - A - k - B   =   m -  C
///       |       |           | |
///       n       l           n l
///
/// NOTE: in case of TArray, this performs many unuse instances of gemv and gemm depend on tensor rank
///
template<
   typename _T,
   class _TensorA, class _TensorB, class _TensorC,
   class _AnnotationA, class _AnnotationB, class _AnnotationC,
   class = typename std::enable_if<
      is_boxtensor<_TensorA>::value &
      is_boxtensor<_TensorB>::value &
      is_boxtensor<_TensorC>::value &
      is_container<_AnnotationA>::value &
      is_container<_AnnotationB>::value &
      is_container<_AnnotationC>::value
   >::type
>
void contract(
   const _T& alpha,
   const _TensorA& A, const _AnnotationA& aA,
   const _TensorB& B, const _AnnotationB& aB,
   const _T& beta,
         _TensorC& C, const _AnnotationC& aC)
{
   // check index A
   auto __sort_indexA = _AnnotationA{aA};
   std::sort(std::begin(__sort_indexA), std::end(__sort_indexA));
   assert(std::unique(std::begin(__sort_indexA), std::end(__sort_indexA)) == std::end(__sort_indexA));

   // check index B
   auto __sort_indexB = _AnnotationB{aB};
   std::sort(std::begin(__sort_indexB), std::end(__sort_indexB));
   assert(std::unique(std::begin(__sort_indexB), std::end(__sort_indexB)) == std::end(__sort_indexB));

   // check index C
   auto __sort_indexC = _AnnotationC{aC};
   std::sort(std::begin(__sort_indexC), std::end(__sort_indexC));
   assert(std::unique(std::begin(__sort_indexC), std::end(__sort_indexC)) == std::end(__sort_indexC));

   typedef btas::varray<size_t> Permutation;

   // permute index A
   Permutation __permute_indexA;
   resize(__permute_indexA, aA.size());

   // permute index B
   Permutation __permute_indexB;
   resize(__permute_indexB, aB.size());

   // permute index C
   Permutation __permute_indexC;
   resize(__permute_indexC, aC.size());

   size_type m = 0;
   size_type n = 0;
   size_type k = 0;

   // row index
   for(auto itrA = std::begin(aA); itrA != std::end(aA); ++itrA)
   {
      if(!std::binary_search(std::begin(__sort_indexB), std::end(__sort_indexB), *itrA))
      {
         __permute_indexA[m] = *itrA;
         __permute_indexC[m] = *itrA;
         ++m;
      }
   }
   // index to be contracted
   for(auto itrA = std::begin(aA); itrA != std::end(aA); ++itrA)
   {
      if( std::binary_search(std::begin(__sort_indexB), std::end(__sort_indexB), *itrA))
      {
         __permute_indexA[m+k] = *itrA;
         __permute_indexB[k]   = *itrA;
         ++k;
      }
   }
   // column index
   for(auto itrB = std::begin(aB); itrB != std::end(aB); ++itrB)
   {
      if(!std::binary_search(std::begin(__sort_indexA), std::end(__sort_indexA), *itrB))
      {
         __permute_indexB[k+n] = *itrB;
         __permute_indexC[m+n] = *itrB;
         ++n;
      }
   }

   // check result index C
   Permutation __sort_permute_indexC(__permute_indexC);
   std::sort(std::begin(__sort_permute_indexC), std::end(__sort_permute_indexC));
   assert(std::equal(std::begin(__sort_permute_indexC), std::end(__sort_permute_indexC), std::begin(__sort_indexC)));

   // permute A
   std::shared_ptr<const _TensorA> __refA;
   if(!std::equal(std::begin(aA), std::end(aA), std::begin(__permute_indexA)))
   {
      __refA = std::make_shared<const _TensorA>();
      permute(A, aA, const_cast<_TensorA&>(*__refA), __permute_indexA);
   }
   else
   {
      __refA.reset(&A, btas::nulldeleter());
   }

   // permute B
   std::shared_ptr<const _TensorB> __refB;
   if(!std::equal(std::begin(aB), std::end(aB), std::begin(__permute_indexB)))
   {
      __refB = std::make_shared<const _TensorB>();
      permute(B, aB, const_cast<_TensorB&>(*__refB), __permute_indexB);
   }
   else
   {
      __refB.reset(&B, btas::nulldeleter());
   }

   bool __C_to_permute = false;

   // to set rank of C
   if(C.empty())
   {
      Permutation __zero_shape;
      resize(__zero_shape, m+n);
      std::fill(std::begin(__zero_shape), std::end(__zero_shape), 0);
      C.resize(__zero_shape);
   }

   // permute C
   std::shared_ptr<_TensorC> __refC;
   if(!std::equal(std::begin(aC), std::end(aC), std::begin(__permute_indexC)))
   {
      __refC = std::make_shared<_TensorC>();
      permute(C, aC, *__refC, __permute_indexC);
      __C_to_permute = true;
   }
   else
   {
      __refC.reset(&C, btas::nulldeleter());
   }

   // call BLAS functions
   if     (rank(A) == k && rank(B) == k)
   {
      assert(false); // dot should be called instead
   }
   else if(k == 0)
   {
      scal(beta, *__refC);
      ger (alpha, *__refA, *__refB, *__refC);
   }
   else if(rank(A) == k)
   {
      gemv(CblasTrans,   alpha, *__refB, *__refA, beta, *__refC);
   }
   else if(rank(B) == k)
   {
      gemv(CblasNoTrans, alpha, *__refA, *__refB, beta, *__refC);
   }
   else
   {
      gemm(CblasNoTrans, CblasNoTrans, alpha, *__refA, *__refB, beta, *__refC);
   }

   // permute back
   if(__C_to_permute)
   {
      permute(*__refC, __permute_indexC, C, aC);
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
      std::is_same<typename _TensorA::value_type, typename _TensorB::value_type>::value &
      std::is_same<typename _TensorA::value_type, typename _TensorC::value_type>::value
   >::type
>
void contract(
   const _T& alpha,
   const _TensorA& A, std::initializer_list<_UA> aA,
   const _TensorB& B, std::initializer_list<_UB> aB,
   const _T& beta,
         _TensorC& C, std::initializer_list<_UC> aC)
{
    contract(alpha,
             A, btas::varray<_UA>(aA),
             B, btas::varray<_UB>(aB),
             beta,
             C, btas::varray<_UC>(aC)
            );
}

} //namespace btas

#endif
