#ifndef __BTAS_GEMV_IMPL_H
#define __BTAS_GEMV_IMPL_H 1

#include <algorithm>
#include <numeric>
#include <iterator>
#include <type_traits>

#include <btas/tensor_traits.h>
#include <btas/types.h>

#include <btas/generic/numeric_type.h>
#include <btas/generic/tensor_iterator_wrapper.h>

#include <btas/generic/scal_impl.h>

namespace btas {


template<bool _Finalize> struct gemv_impl { };

template<> struct gemv_impl<true>
{
#if 0
   /// GEMV implementation
   template<typename _T, class _IteratorA, class _IteratorX, class _IteratorY>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorA itrA,
      const unsigned long& LDA,
            _IteratorX itrX,
      const typename std::iterator_traits<_IteratorX>::difference_type& incX,
      const _T& beta,
            _IteratorY itrY,
      const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      if (beta != NumericType<_T>::one())
      {
         if (transA == CblasNoTrans)
            scal (Msize, beta, itrY, incY);
         else
            scal (Nsize, beta, itrY, incY);
      }

      // A:NoTrans RowMajor
      if      (transA == CblasNoTrans && order == CblasRowMajor)
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Msize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrA, ++itrX)
            {
               (*itrY) += alpha * (*itrA) * (*itrX);
            }
         }
      }
      // A:Trans RowMajor
      else if (transA == CblasTrans && order == CblasRowMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Msize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrA, ++itrY)
            {
               (*itrY) += alpha * (*itrA) * (*itrX);
            }
         }
      }
      // A:ConjTrans RowMajor
      else if (transA == CblasConjTrans && order == CblasRowMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Msize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrA, ++itrY)
            {
               (*itrY) += alpha * impl::conj(*itrA) * (*itrX);
            }
         }
      }
      // A:NoTrans ColMajor
      else if (transA == CblasNoTrans && order == CblasColMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Nsize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Msize; ++j, ++itrA, ++itrY)
            {
               (*itrY) += alpha * (*itrA) * (*itrX);
            }
         }
      }
      // A:Trans ColMajor
      else if (transA == CblasTrans && order == CblasColMajor)
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Nsize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Msize; ++j, ++itrA, ++itrX)
            {
               (*itrY) += alpha * (*itrA) * (*itrX);
            }
         }
      }
      // A:ConjTrans ColMajor
      else if (transA == CblasConjTrans && order == CblasColMajor)
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Nsize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Msize; ++j, ++itrA, ++itrX)
            {
               (*itrY) += alpha * impl::conj(*itrA) * (*itrX);
            }
         }
      }
   }
#endif

#ifdef _HAS_CBLAS

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, float>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const float* itrA,
      const unsigned long& LDA,
      const float* itrX,
      const typename std::iterator_traits<float*>::difference_type& incX,
      const _T& beta,
            float* itrY,
      const typename std::iterator_traits<float*>::difference_type& incY)
   {
      cblas_sgemv(order, transA, Msize, Nsize, alpha, itrA, LDA, itrX, incX, beta, itrY, incY);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, double>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const double* itrA,
      const unsigned long& LDA,
      const double* itrX,
      const typename std::iterator_traits<double*>::difference_type& incX,
      const _T& beta,
            double* itrY,
      const typename std::iterator_traits<double*>::difference_type& incY)
   {
      cblas_dgemv(order, transA, Msize, Nsize, alpha, itrA, LDA, itrX, incX, beta, itrY, incY);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, std::complex<float>>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const std::complex<float>* itrA,
      const unsigned long& LDA,
      const std::complex<float>* itrX,
      const typename std::iterator_traits<std::complex<float>*>::difference_type& incX,
      const _T& beta,
            std::complex<float>* itrY,
      const typename std::iterator_traits<std::complex<float>*>::difference_type& incY)
   {
      const std::complex<float> alphac(std::move(alpha));
      const std::complex<float> betac (std::move(beta));
      cblas_cgemv(order, transA, Msize, Nsize, &alphac, itrA, LDA, itrX, incX, &betac, itrY, incY);
   }

   template <typename _T, class = typename std::enable_if<std::is_convertible<_T, std::complex<double>>::value>::type>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
      const std::complex<double>* itrA,
      const unsigned long& LDA,
      const std::complex<double>* itrX,
      const typename std::iterator_traits<std::complex<double>*>::difference_type& incX,
      const _T& beta,
            std::complex<double>* itrY,
      const typename std::iterator_traits<std::complex<double>*>::difference_type& incY)
   {
      const std::complex<double> alphac(std::move(alpha));
      const std::complex<double> betac (std::move(beta));
      cblas_zgemv(order, transA, Msize, Nsize, &alphac, itrA, LDA, itrX, incX, &betac, itrY, incY);
   }

#endif // _HAS_CBLAS

};

template<> struct gemv_impl<false>
{
   /// GEMV implementation
   template<typename _T, class _IteratorA, class _IteratorX, class _IteratorY>
   static void call (
      const CBLAS_ORDER& order,
      const CBLAS_TRANSPOSE& transA,
      const unsigned long& Msize,
      const unsigned long& Nsize,
      const _T& alpha,
            _IteratorA itrA,
      const unsigned long& LDA,
            _IteratorX itrX,
      const typename std::iterator_traits<_IteratorX>::difference_type& incX,
      const _T& beta,
            _IteratorY itrY,
      const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      // A:NoTrans RowMajor
      if      (transA == CblasNoTrans && order == CblasRowMajor)
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Msize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrA, ++itrX)
            {
               gemv(order, transA, alpha, *itrA, *itrX, beta, *itrY);
            }
         }
      }
      // A:Trans RowMajor
      else if (transA != CblasNoTrans && order == CblasRowMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Msize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Nsize; ++j, ++itrA, ++itrY)
            {
               gemv(order, transA, alpha, *itrA, *itrX, beta, *itrY);
            }
         }
      }
      // A:NoTrans ColMajor
      else if (transA == CblasNoTrans && order == CblasColMajor)
      {
         auto itrY_save = itrY;
         for (size_type i = 0; i < Nsize; ++i, ++itrX)
         {
            itrY = itrY_save;
            for (size_type j = 0; j < Msize; ++j, ++itrA, ++itrY)
            {
               gemv(order, transA, alpha, *itrA, *itrX, beta, *itrY);
            }
         }
      }
      // A:Trans ColMajor
      else if (transA != CblasNoTrans && order == CblasColMajor)
      {
         auto itrX_save = itrX;
         for (size_type i = 0; i < Nsize; ++i, ++itrY)
         {
            itrX = itrX_save;
            for (size_type j = 0; j < Msize; ++j, ++itrA, ++itrX)
            {
               gemv(order, transA, alpha, *itrA, *itrX, beta, *itrY);
            }
         }
      }
   }

};

//  ================================================================================================

/// Generic implementation of BLAS GEMV in terms of C++ iterator
template<typename _T, class _IteratorA, class _IteratorX, class _IteratorY>
void gemv (
   const CBLAS_ORDER& order,
   const CBLAS_TRANSPOSE& transA,
   const unsigned long& Msize,
   const unsigned long& Nsize,
   const _T& alpha,
         _IteratorA itrA,
   const unsigned long& LDA,
         _IteratorX itrX,
   const typename std::iterator_traits<_IteratorX>::difference_type& incX,
   const _T& beta,
         _IteratorY itrY,
   const typename std::iterator_traits<_IteratorY>::difference_type& incY)
{
   typedef std::iterator_traits<_IteratorA> __traits_A;
   typedef std::iterator_traits<_IteratorX> __traits_X;
   typedef std::iterator_traits<_IteratorY> __traits_Y;

   typedef typename __traits_A::value_type value_type;

   static_assert(std::is_same<value_type, typename __traits_X::value_type>::value, "value type of X must be the same as that of A");
   static_assert(std::is_same<value_type, typename __traits_Y::value_type>::value, "value type of Y must be the same as that of A");

   static_assert(std::is_same<typename __traits_A::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator A must be a random access iterator");

   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator X must be a random access iterator");

   static_assert(std::is_same<typename __traits_Y::iterator_category, std::random_access_iterator_tag>::value,
                 "iterator Y must be a random access iterator");

   typename _IteratorA::pointer A = &(*itrA);
   typename _IteratorX::pointer X = &(*itrX);
   typename _IteratorY::pointer Y = &(*itrY);
   gemv_impl<std::is_convertible<_T, value_type>::value>::call(order, transA, Msize, Nsize, alpha, A, LDA, X, incX, beta, Y, incY);
   //gemv_impl<std::is_convertible<_T, value_type>::value>::call(order, transA, Msize, Nsize, alpha, itrA, LDA, itrX, incX, beta, itrY, incY);
}

//  ================================================================================================

/// Generic interface of BLAS-GEMV
/// \param order storage order of tensor in matrix view (CblasRowMajor, CblasColMajor)
/// \param transA transpose directive for tensor A (CblasNoTrans, CblasTrans, CblasConjTrans)
/// \param alpha scalar value to be multiplied to A * X
/// \param A input tensor
/// \param X input tensor
/// \param beta scalar value to be multiplied to Y
/// \param Y output tensor which can be empty tensor but needs to have rank info (= size of shape).
/// Iterator is assumed to be consecutive (or, random_access_iterator) , thus e.g. iterator to map doesn't work.
template<
   typename _T,
   class _TensorA, class _TensorX, class _TensorY,
   class = typename std::enable_if<
      is_boxtensor<_TensorA>::value &
      is_boxtensor<_TensorX>::value &
      is_boxtensor<_TensorY>::value
   >::type
>
void gemv (
   const CBLAS_TRANSPOSE& transA,
   const _T& alpha,
   const _TensorA& A,
   const _TensorX& X,
   const _T& beta,
         _TensorY& Y)
{
    static_assert(boxtensor_storage_order<_TensorA>::value == boxtensor_storage_order<_TensorY>::value &&
                  boxtensor_storage_order<_TensorX>::value == boxtensor_storage_order<_TensorY>::value,
                  "btas::gemm does not support mixed storage order");
    const CBLAS_ORDER order = boxtensor_storage_order<_TensorY>::value == boxtensor_storage_order<_TensorY>::row_major ?
                              CblasRowMajor : CblasColMajor;

   if (A.empty() || X.empty())
   {
      scal(beta, Y);
      return;
   }

   assert(not ((transA == CblasConjTrans ) && std::is_fundamental<typename _TensorA::value_type>::value));

   // get contraction rank
   const size_type rankX = rank(X);
   const size_type rankY = rank(Y);

   // get shapes
   const typename _TensorA::range_type::extent_type& extentA = extent(A);
   const typename _TensorX::range_type::extent_type& extentX = extent(X);
         typename _TensorY::range_type::extent_type  extentY = extent(Y); // if Y is empty, this gives { 0,...,0 }

   size_type Msize = 0; // Rows count of Y
   size_type Nsize = 0; // Cols count of Y

   size_type LDA = 0; // Leading dims of A

   // to minimize forks by if?
   if (transA == CblasNoTrans)
   {
      Msize = std::accumulate(std::begin(extentA), std::begin(extentA)+rankY, 1ul, std::multiplies<size_type>());
      Nsize = std::accumulate(std::begin(extentA)+rankY, std::end(extentA),   1ul, std::multiplies<size_type>());

      for (size_type i = 0; i < rankY; ++i) extentY[i] = extentA[i];

      assert(std::equal(std::begin(extentA)+rankY, std::end(extentA), std::begin(extentX)));
   }
   else
   {
      Msize = std::accumulate(std::begin(extentA), std::begin(extentA)+rankX, 1ul, std::multiplies<size_type>());
      Nsize = std::accumulate(std::begin(extentA)+rankX, std::end(extentA),   1ul, std::multiplies<size_type>());

      for (size_type i = 0; i < rankY; ++i) extentY[i]   = extentA[i+rankX];

      assert(std::equal(std::begin(extentA), std::begin(extentA)+rankX, std::begin(extentX)));
   }

   LDA = std::accumulate(std::begin(extentA)+rankY, std::end(extentA),   1ul, std::multiplies<size_type>());

   //if(order == CblasRowMajor)
   //{
   //   LDA = Nsize;
   //}
   //else
   //{
   //   LDA = Msize;
   //}

   // resize / scale
   if (Y.empty())
   {
     typedef typename _TensorY::value_type value_type;
      Y.resize(extentY);
      NumericType<value_type>::fill(std::begin(Y), std::end(Y), NumericType<value_type>::zero());
   }
   else
   {
      assert(std::equal(std::begin(extentY), std::end(extentY), std::begin(extent(Y))));
   }

   auto itrA = std::begin(A);
   auto itrX = std::begin(X);
   auto itrY = std::begin(Y);

   gemv (order, transA, Msize, Nsize, alpha, itrA, LDA, itrX, 1, beta, itrY, 1);
}

} // namespace btas

#endif // __BTAS_GEMV_IMPL_H
