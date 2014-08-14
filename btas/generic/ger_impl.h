#ifndef __BTAS_GER_IMPL_H
#define __BTAS_GER_IMPL_H 1

#include <algorithm>
#include <numeric>
#include <iterator>
#include <type_traits>

#include <btas/tensor_traits.h>
#include <btas/types.h>

#include <btas/generic/numeric_type.h>
#include <btas/generic/tensor_iterator_wrapper.h>

namespace btas {

template<bool _Finalize> struct ger_impl { };

template<> struct ger_impl<true>
{
   /// Performs GER operation
   template<typename _T, class _IteratorX, class _IteratorY, class _IteratorA>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX,
      const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY,
      const typename std::iterator_traits<_IteratorY>::difference_type& incY,
            _IteratorA itrA,
      const unsigned long& LDA)
   {
      // RowMajor
      if (order == CblasRowMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Msize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrY, ++itrA)
            {
               (*itrA) += alpha * (*itrX) * (*itrY);
            }
         }
      }
      // A: ColMajor
      else
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Nsize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Msize; ++j, ++itrX, ++itrA)
            {
               (*itrA) += alpha * (*itrX) * (*itrY);
            }
         }
      }
   }

#ifdef _HAS_CBLAS

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, float>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const float* itrX,
      const typename std::iterator_traits<float*>::difference_type& incX,
      const float* itrY,
      const typename std::iterator_traits<float*>::difference_type& incY,
            float* itrA,
      const unsigned long& LDA)
   {
      cblas_sger(order, Msize, Nsize, alpha, itrX, incX, itrY, incY, itrA, LDA);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, double>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const double* itrX,
      const typename std::iterator_traits<double*>::difference_type& incX,
      const double* itrY,
      const typename std::iterator_traits<double*>::difference_type& incY,
            double* itrA,
      const unsigned long& LDA)
   {
      cblas_dger(order, Msize, Nsize, alpha, itrX, incX, itrY, incY, itrA, LDA);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, std::complex<float>>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const std::complex<float>* itrX,
      const typename std::iterator_traits<std::complex<float>*>::difference_type& incX,
      const std::complex<float>* itrY,
      const typename std::iterator_traits<std::complex<float>*>::difference_type& incY,
            std::complex<float>* itrA,
      const unsigned long& LDA)
   {
      // FIXME: implement cgerc and cgeru separately.
      const std::complex<float> alphac(std::move(alpha));
      cblas_cgeru(order, Msize, Nsize, &alphac, itrX, incX, itrY, incY, itrA, LDA);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, std::complex<double>>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const std::complex<double>* itrX,
      const typename std::iterator_traits<std::complex<double>*>::difference_type& incX,
      const std::complex<double>* itrY,
      const typename std::iterator_traits<std::complex<double>*>::difference_type& incY,
            std::complex<double>* itrA,
      const unsigned long& LDA)
   {
      // FIXME: implement zgerc and zgeru separately.
      const std::complex<double> alphac(std::move(alpha));
      cblas_zgeru(order, Msize, Nsize, &alphac, itrX, incX, itrY, incY, itrA, LDA);
   }

#endif // _HAS_CBLAS

};

template<> struct ger_impl<false>
{
   /// GER implementation
   template<typename _T, class _IteratorX, class _IteratorY, class _IteratorA>
   static void call (
      const CBLAS_ORDER& order,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorX itrX,
      const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY,
      const typename std::iterator_traits<_IteratorY>::difference_type& incY,
            _IteratorA itrA,
      const unsigned long& LDA)
   {
      // RowMajor
      if (order == CblasRowMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Msize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrY, ++itrA)
            {
               ger(order, alpha, *itrX, *itrY, *itrA);
            }
         }
      }
      // A: ColMajor
      else
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Nsize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Msize; ++j, ++itrX, ++itrA)
            {
               ger(order, alpha, *itrX, *itrY, *itrA);
            }
         }
      }
   }

};

//  ================================================================================================

/// Generic implementation of BLAS GER in terms of C++ iterator
template<typename _T, class _IteratorX, class _IteratorY, class _IteratorA>
void ger (
   const CBLAS_ORDER& order,
   const unsigned long& Msize,
   const unsigned long& Nsize,
   const _T& alpha,
         _IteratorX itrX,
   const typename std::iterator_traits<_IteratorX>::difference_type& incX,
         _IteratorY itrY,
   const typename std::iterator_traits<_IteratorY>::difference_type& incY,
         _IteratorA itrA,
   const unsigned long& LDA)
{
   typedef std::iterator_traits<_IteratorX> __traits_X;
   typedef std::iterator_traits<_IteratorY> __traits_Y;
   typedef std::iterator_traits<_IteratorA> __traits_A;

   typedef typename __traits_A::value_type value_type;

   static_assert(std::is_same<value_type, typename __traits_X::value_type>::value, "value type of X must be the same as that of A");
   static_assert(std::is_same<value_type, typename __traits_Y::value_type>::value, "value type of Y must be the same as that of A");

   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator X must be a random access iterator");

   static_assert(std::is_same<typename __traits_Y::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator Y must be a random access iterator");

   static_assert(std::is_same<typename __traits_A::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator A must be a random access iterator");

   ger_impl<std::is_convertible<_T, value_type>::value>::call(order, Msize, Nsize, alpha, itrX, incX, itrY, incY, itrA, LDA);
}

//  ================================================================================================

/// Generic implementation of operation (generalization of BLAS GER operation)
/// \param alpha scalar value to be multiplied to A * X * Y
/// \param X input tensor
/// \param Y input tensor
/// \param A output tensor which can be empty tensor but needs to have rank info (= size of shape).
/// Iterator is assumed to be consecutive (or, random_access_iterator) , thus e.g. iterator to map doesn't work.
template<
   typename _T,
   class _TensorX, class _TensorY, class _TensorA,
   class = typename std::enable_if<
      is_boxtensor<_TensorX>::value &
      is_boxtensor<_TensorY>::value &
      is_boxtensor<_TensorA>::value
   >::type
>
void ger (
   const _T& alpha,
   const _TensorX& X,
   const _TensorY& Y,
         _TensorA& A)
{
    static_assert(boxtensor_storage_order<_TensorX>::value == boxtensor_storage_order<_TensorA>::value &&
                  boxtensor_storage_order<_TensorY>::value == boxtensor_storage_order<_TensorA>::value,
                  "btas::ger does not support mixed storage order");
    const CBLAS_ORDER order = boxtensor_storage_order<_TensorA>::value == boxtensor_storage_order<_TensorA>::row_major ?
                              CblasRowMajor : CblasColMajor;

   if (X.empty() || Y.empty())
   {
      return;
   }

   // get contraction rank
   const size_type rankX = rank(X);
   const size_type rankY = rank(Y);

   // get shapes
   const typename _TensorX::range_type::extent_type& extentX = extent(X);
   const typename _TensorY::range_type::extent_type& extentY = extent(Y);
         typename _TensorA::range_type::extent_type  extentA = extent(A);

   size_type Msize = std::accumulate(std::begin(extentX), std::end(extentX), 1ul, std::multiplies<size_type>());
   size_type Nsize = std::accumulate(std::begin(extentY), std::end(extentY), 1ul, std::multiplies<size_type>());
   size_type LDA   = (order == CblasRowMajor) ? Nsize : Msize;

   std::copy_n(std::begin(extentX), rankX, std::begin(extentA));
   std::copy_n(std::begin(extentY), rankY, std::begin(extentA)+rankX);

   // resize / scale
   if (A.empty())
   {
     typedef typename _TensorA::value_type value_type;
     A.resize(extentA);
     NumericType<value_type>::fill(std::begin(A), std::end(A), NumericType<value_type>::zero());
   }
   else
   {
      assert(std::equal(std::begin(extentA), std::end(extentA), std::begin(extent(A))));
   }

   auto itrX = std::begin(X);
   auto itrY = std::begin(Y);
   auto itrA = std::begin(A);

   ger (order, Msize, Nsize, alpha, itrX, 1, itrY, 1, itrA, LDA);
}

} // namespace btas

#endif // __BTAS_GER_IMPL_H
