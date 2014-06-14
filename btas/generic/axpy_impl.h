#ifndef __BTAS_AXPY_IMPL_H
#define __BTAS_AXPY_IMPL_H 1

#include <algorithm>
#include <iterator>
#include <type_traits>

#include <btas/tensor.h>
#include <btas/tensor_traits.h>
#include <btas/types.h>

#include <btas/generic/numeric_type.h>
#include <btas/generic/tensor_iterator_wrapper.h>

namespace btas {

//  ================================================================================================

/// Call BLAS depending on type of Tensor class
template<bool _Finalize> struct axpy_impl { };

/// Case that alpha is trivially multipliable to elements
template<> struct axpy_impl<true>
{
   template<typename _T, class _IteratorX, class _IteratorY>
   static void call (
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      for (unsigned long i = 0; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         (*itrY) += alpha * (*itrX);
      }
   }

#ifdef _HAS_CBLAS
   static void call (
      const unsigned long& Nsize,
      const float& alpha,
      const float* itrX, const typename std::iterator_traits<float*>::difference_type& incX,
            float* itrY, const typename std::iterator_traits<float*>::difference_type& incY)
   {
      cblas_saxpy(Nsize, alpha, itrX, incX, itrY, incY);
   }

   static void call (
      const unsigned long& Nsize,
      const double& alpha,
      const double* itrX, const typename std::iterator_traits<double*>::difference_type& incX,
            double* itrY, const typename std::iterator_traits<double*>::difference_type& incY)
   {
      cblas_daxpy(Nsize, alpha, itrX, incX, itrY, incY);
   }

   static void call (
      const unsigned long& Nsize,
      const std::complex<float>& alpha,
      const std::complex<float>* itrX, const typename std::iterator_traits<std::complex<float>*>::difference_type& incX,
            std::complex<float>* itrY, const typename std::iterator_traits<std::complex<float>*>::difference_type& incY)
   {
      cblas_caxpy(Nsize, &alpha, itrX, incX, itrY, incY);
   }

   static void call (
      const unsigned long& Nsize,
      const std::complex<double>& alpha,
      const std::complex<double>* itrX, const typename std::iterator_traits<std::complex<double>*>::difference_type& incX,
            std::complex<double>* itrY, const typename std::iterator_traits<std::complex<double>*>::difference_type& incY)
   {
      cblas_zaxpy(Nsize, &alpha, itrX, incX, itrY, incY);
   }
#endif
};

/// Case that alpha is multiplied recursively by AXPY
/// Note that incX and incY are disabled for recursive call
template<> struct axpy_impl<false>
{
   template<typename _T, class _IteratorX, class _IteratorY>
   static void call (
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      for (unsigned long i = 0; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         axpy(alpha, *itrX, *itrY);
      }
   }
};

//  ================================================================================================

/// Generic implementation of BLAS AXPY in terms of C++ iterator
template<typename _T, class _IteratorX, class _IteratorY>
void axpy (
   const unsigned long& Nsize,
   const _T& alpha,
         _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
         _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
{
   typedef std::iterator_traits<_IteratorX> __traits_X;
   typedef std::iterator_traits<_IteratorY> __traits_Y;

   static_assert(std::is_same<typename __traits_X::value_type, typename __traits_Y::value_type>::value, "value type of Y must be the same as that of X");
   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value, "iterator X must be a random access iterator");
   static_assert(std::is_same<typename __traits_Y::iterator_category, std::random_access_iterator_tag>::value, "iterator Y must be a random access iterator");

   const bool match = std::is_convertible<_T, typename __traits_X::value_type>::value;
   axpy_impl<std::is_convertible<_T, typename __traits_X::value_type>::value>
          ::call(Nsize, !match ? alpha : static_cast<typename __traits_X::value_type>(alpha), itrX, incX, itrY, incY);
}

//  ================================================================================================

/// Convenient wrapper to call BLAS AXPY from tensor objects
template<
   typename _T,
   class _TensorX, class _TensorY,
   class = typename std::enable_if<
      is_boxtensor<_TensorX>::value &
      is_boxtensor<_TensorY>::value
   >::type
>
void axpy (
   const _T& alpha,
   const _TensorX& X,
         _TensorY& Y)
{
   typedef typename _TensorX::value_type value_type;
   static_assert(std::is_same<value_type, typename _TensorY::value_type>::value, "value type of Y must be the same as that of X");

   if (X.empty())
   {
      Y.clear();
      return;
   }

   if (Y.empty())
   {
      Y.resize(btas::extent(X));
      NumericType<value_type>::fill(std::begin(Y), std::end(Y), NumericType<value_type>::zero());
   }
   else
   {
      assert( range(X) == range(Y) );
   }

   auto itrX = std::begin(X);
   auto itrY = std::begin(Y);

   axpy (X.size(), alpha, itrX, 1, itrY, 1);
}

} // namespace btas

#endif // __BTAS_AXPY_IMPL_H
