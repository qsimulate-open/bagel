#ifndef __BTAS_DOT_IMPL_H
#define __BTAS_DOT_IMPL_H 1

#include <algorithm>
#include <iterator>
#include <type_traits>

#include <btas/tensor_traits.h>
#include <btas/types.h>

#include <btas/generic/numeric_type.h>
#include <btas/generic/tensor_iterator_wrapper.h>

namespace btas {

template<typename _T>
struct __dot_result_type {
   typedef typename __dot_result_type<typename _T::value_type>::type type;
};

template<> struct __dot_result_type<float> { typedef float type; };

template<> struct __dot_result_type<double> { typedef double type; };

template<> struct __dot_result_type<std::complex<float>> { typedef std::complex<float> type; };

template<> struct __dot_result_type<std::complex<double>> { typedef std::complex<double> type; };

//  ================================================================================================

/// Dot with conjugation for general case
template<typename _T>
struct dotc_impl
{
   typedef typename __dot_result_type<_T>::type return_type;

   template<class _IteratorX, class _IteratorY>
   static return_type call (
      const unsigned long& Nsize,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      return_type val = dotc(*itrX, *itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += dotc(*itrX, *itrY);
      }
      return val;
   }
};

/// Dot without conjugation for general case
template<typename _T>
struct dotu_impl
{
   typedef typename __dot_result_type<_T>::type return_type;

   template<class _IteratorX, class _IteratorY>
   static return_type call (
      const unsigned long& Nsize,
            _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
            _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
   {
      // redirect to dotc, otherwise must be specialized
      return dotc_impl<_T>::call(Nsize, itrX, incX, itrY, incY);
   }
};

template<>
struct dotc_impl<float>
{
   typedef float return_type;

   static return_type call (
      const unsigned long& Nsize,
      const float* itrX, const typename std::iterator_traits<float*>::difference_type& incX,
      const float* itrY, const typename std::iterator_traits<float*>::difference_type& incY)
   {
#ifdef _HAS_CBLAS
      return cblas_sdot(Nsize, itrX, incX, itrY, incY);
#else
      return_type val = (*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += (*itrX) * (*itrY);
      }
      return val;
#endif
   }
};

template<>
struct dotc_impl<double>
{
   typedef double return_type;

   static return_type call (
      const unsigned long& Nsize,
      const double* itrX, const typename std::iterator_traits<double*>::difference_type& incX,
      const double* itrY, const typename std::iterator_traits<double*>::difference_type& incY)
   {
#ifdef _HAS_CBLAS
      return cblas_ddot(Nsize, itrX, incX, itrY, incY);
#else
      return_type val = (*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += (*itrX) * (*itrY);
      }
      return val;
#endif
   }
};

template<>
struct dotc_impl<std::complex<float>>
{
   typedef std::complex<float> return_type;

   static return_type call (
      const unsigned long& Nsize,
      const std::complex<float>* itrX, const typename std::iterator_traits<std::complex<float>*>::difference_type& incX,
      const std::complex<float>* itrY, const typename std::iterator_traits<std::complex<float>*>::difference_type& incY)
   {
      return_type val;
#ifdef _HAS_CBLAS
      cblas_cdotc_sub(Nsize, itrX, incX, itrY, incY, &val);
#else
      val = std::conj(*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += std::conj(*itrX) * (*itrY);
      }
#endif
      return val;
   }
};

template<>
struct dotu_impl<std::complex<float>>
{
   typedef std::complex<float> return_type;

   static return_type call (
      const unsigned long& Nsize,
      const std::complex<float>* itrX, const typename std::iterator_traits<std::complex<float>*>::difference_type& incX,
      const std::complex<float>* itrY, const typename std::iterator_traits<std::complex<float>*>::difference_type& incY)
   {
      return_type val;
#ifdef _HAS_CBLAS
      cblas_cdotu_sub(Nsize, itrX, incX, itrY, incY, &val);
#else
      val = (*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += *itrX * (*itrY);
      }
#endif
      return val;
   }
};

template<>
struct dotc_impl<std::complex<double>>
{
   typedef std::complex<double> return_type;

   static return_type call (
      const unsigned long& Nsize,
      const std::complex<double>* itrX, const typename std::iterator_traits<std::complex<double>*>::difference_type& incX,
      const std::complex<double>* itrY, const typename std::iterator_traits<std::complex<double>*>::difference_type& incY)
   {
      return_type val;
#ifdef _HAS_CBLAS
      cblas_zdotc_sub(Nsize, itrX, incX, itrY, incY, &val);
#else
      val = std::conj(*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += std::conj(*itrX) * (*itrY);
      }
#endif
      return val;
   }
};

template<>
struct dotu_impl<std::complex<double>>
{
   typedef std::complex<double> return_type;

   static return_type call (
      const unsigned long& Nsize,
      const std::complex<double>* itrX, const typename std::iterator_traits<std::complex<double>*>::difference_type& incX,
      const std::complex<double>* itrY, const typename std::iterator_traits<std::complex<double>*>::difference_type& incY)
   {
      return_type val;
#ifdef _HAS_CBLAS
      cblas_zdotu_sub(Nsize, itrX, incX, itrY, incY, &val);
#else
      val = (*itrX) * (*itrY);
      itrX += incX;
      itrY += incY;
      for (unsigned long i = 1; i < Nsize; ++i, itrX += incX, itrY += incY)
      {
         val += *itrX * (*itrY);
      }
#endif
      return val;
   }
};

//  ================================================================================================

/// Generic implementation of BLAS DOT in terms of C++ iterator
template<class _IteratorX, class _IteratorY>
typename __dot_result_type<typename std::iterator_traits<_IteratorX>::value_type>::type
dotc (
   const unsigned long& Nsize,
         _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
         _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
{
   typedef std::iterator_traits<_IteratorX> __traits_X;
   typedef std::iterator_traits<_IteratorY> __traits_Y;

   static_assert(std::is_same<typename __traits_X::value_type, typename __traits_Y::value_type>::value, "value type of Y must be the same as that of X");
   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value, "iterator X must be a random access iterator");
   static_assert(std::is_same<typename __traits_Y::iterator_category, std::random_access_iterator_tag>::value, "iterator Y must be a random access iterator");

   return dotc_impl<typename __traits_X::value_type>::call(Nsize, itrX, incX, itrY, incY);
}

/// Generic implementation of BLAS DOT in terms of C++ iterator
template<class _IteratorX, class _IteratorY>
typename __dot_result_type<typename std::iterator_traits<_IteratorX>::value_type>::type
dotu (
   const unsigned long& Nsize,
         _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
         _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
{
   typedef std::iterator_traits<_IteratorX> __traits_X;
   typedef std::iterator_traits<_IteratorY> __traits_Y;

   static_assert(std::is_same<typename __traits_X::value_type, typename __traits_Y::value_type>::value, "value type of Y must be the same as that of X");
   static_assert(std::is_same<typename __traits_X::iterator_category, std::random_access_iterator_tag>::value, "iterator X must be a random access iterator");
   static_assert(std::is_same<typename __traits_Y::iterator_category, std::random_access_iterator_tag>::value, "iterator Y must be a random access iterator");

   return dotu_impl<typename __traits_X::value_type>::call(Nsize, itrX, incX, itrY, incY);
}

/// Generic implementation of BLAS DOT in terms of C++ iterator
template<class _IteratorX, class _IteratorY>
typename __dot_result_type<typename std::iterator_traits<_IteratorX>::value_type>::type
dot (
   const unsigned long& Nsize,
         _IteratorX itrX, const typename std::iterator_traits<_IteratorX>::difference_type& incX,
         _IteratorY itrY, const typename std::iterator_traits<_IteratorY>::difference_type& incY)
{
   return dotc(Nsize, itrX, incX, itrY, incY);
}

//  ================================================================================================

/// Convenient wrapper to call BLAS DOT-C from tensor objects
template<
   class _TensorX,
   class _TensorY,
   class = typename std::enable_if<
      is_tensor<_TensorX>::value &
      is_tensor<_TensorY>::value
   >::type
>
typename __dot_result_type<typename _TensorX::value_type>::type
dotc (const _TensorX& X, const _TensorY& Y)
{
   typedef typename _TensorX::value_type value_type;
   static_assert(std::is_same<value_type, typename _TensorY::value_type>::value, "value type of Y must be the same as that of X");

   if (X.empty() || Y.empty())
   {
      return 0;
   }

   auto itrX = tbegin(X);
   auto itrY = tbegin(Y);

   return dotc(X.size(), itrX, 1, itrY, 1);
}

/// Convenient wrapper to call BLAS DOT-U from tensor objects
template<
   class _TensorX,
   class _TensorY,
   class = typename std::enable_if<
      is_tensor<_TensorX>::value &
      is_tensor<_TensorY>::value
   >::type
>
typename __dot_result_type<typename _TensorX::value_type>::type
dotu (const _TensorX& X, const _TensorY& Y)
{
   typedef typename _TensorX::value_type value_type;
   static_assert(std::is_same<value_type, typename _TensorY::value_type>::value, "value type of Y must be the same as that of X");

   if (X.empty() || Y.empty())
   {
      return 0;
   }

   auto itrX = tbegin(X);
   auto itrY = tbegin(Y);

   return dotu(X.size(), itrX, 1, itrY, 1);
}

/// Convenient wrapper to call BLAS DOT from tensor objects
template<
   class _TensorX,
   class _TensorY,
   class = typename std::enable_if<
      is_tensor<_TensorX>::value &
      is_tensor<_TensorY>::value
   >::type
>
typename __dot_result_type<typename _TensorX::value_type>::type
dot (const _TensorX& X, const _TensorY& Y)
{
   return dotc(X, Y);
}

} // namespace btas

#endif // __BTAS_DOT_IMPL_H
