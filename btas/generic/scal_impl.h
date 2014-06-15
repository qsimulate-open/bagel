#ifndef __BTAS_SCAL_IMPL_H
#define __BTAS_SCAL_IMPL_H 1

#include <algorithm>
#include <iterator>
#include <type_traits>

#include <btas/tensor_traits.h>
#include <btas/types.h>

#include <btas/generic/numeric_type.h>
#include <btas/generic/tensor_iterator_wrapper.h>

namespace btas {

//  ================================================================================================

/// Call BLAS depending on type of Tensor class
template<bool _Finalize> struct scal_impl { };

/// Case that alpha is trivially multipliable to elements
template<> struct scal_impl<true>
{
   template<typename _T, class _IteratorX>
   static void call (
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX)
   {
      for (unsigned long i = 0; i < Nsize; ++i, itrX += incX)
      {
         (*itrX) *= alpha;
      }
   }

#ifdef _HAS_CBLAS
   static void call (
      const unsigned long& Nsize,
      const float& alpha,
            float* itrX, const typename std::iterator_traits<float*>::difference_type& incX)
   {
      cblas_sscal(Nsize, alpha, itrX, incX);
   }

   static void call (
      const unsigned long& Nsize,
      const double& alpha,
            double* itrX, const typename std::iterator_traits<double*>::difference_type& incX)
   {
      cblas_dscal(Nsize, alpha, itrX, incX);
   }

   static void call (
      const unsigned long& Nsize,
      const std::complex<float>& alpha,
            std::complex<float>* itrX, const typename std::iterator_traits<std::complex<float>*>::difference_type& incX)
   {
      cblas_cscal(Nsize, &alpha, itrX, incX);
   }

   static void call (
      const unsigned long& Nsize,
      const std::complex<double>& alpha,
            std::complex<double>* itrX, const typename std::iterator_traits<std::complex<double>*>::difference_type& incX)
   {
      cblas_zscal(Nsize, &alpha, itrX, incX);
   }
#endif
};

/// Case that alpha is multiplied recursively by SCAL
/// Note that incX is disabled for recursive call
template<> struct scal_impl<false>
{
   template<typename _T, class _IteratorX>
   static void call (
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX)
   {
      for (unsigned long i = 0; i < Nsize; ++i, itrX += incX)
      {
         scal(alpha, *itrX);
      }
   }
};

//  ================================================================================================

/// Generic implementation of BLAS SCAL in terms of C++ iterator
template<typename _T, class _IteratorX>
void scal (
   const unsigned long& Nsize,
   const _T& alpha,
         _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX)
{
   typedef std::iterator_traits<_IteratorX> __traits_X;

   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value, "iterator X must be a random access iterator");

   typedef typename __traits_X::value_type __value_X;
   typedef typename std::conditional<std::is_convertible<_T, __value_X>::value, __value_X, _T>::type __alpha;
   scal_impl<std::is_convertible<_T, __value_X>::value>::call(Nsize, static_cast<__alpha>(alpha), itrX, incX);
}

//  ================================================================================================

/// Convenient wrapper to call BLAS SCAL from tensor objects
template<
   typename _T,
   class _TensorX,
   class = typename std::enable_if<
      is_boxtensor<_TensorX>::value
   >::type
>
void scal (
   const _T& alpha,
         _TensorX& X)
{
   if (X.empty())
   {
      return;
   }

   auto itrX = std::begin(X);

   scal (X.size(), alpha, itrX, 1);
}

} // namespace btas

#endif // __BTAS_SCAL_IMPL_H
