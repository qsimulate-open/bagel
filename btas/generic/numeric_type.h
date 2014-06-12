#ifndef __BTAS_NUMERIC_TYPE_H
#define __BTAS_NUMERIC_TYPE_H 1

#include <complex>
#include <algorithm>
#include <iterator>
#include <type_traits>

namespace btas {

/// Numeric value functions
template<typename _T> struct NumericType
{
   static int zero () { return 0; }
   static int one  () { return 1; }

   template<class _Iterator>
   static void fill(_Iterator, _Iterator, int) { }

   template<class _Iterator>
   static void scal(_Iterator, _Iterator, int) { }
};

//
// Specialization for each numeric value type
//

/// Single precision real number
template <> struct NumericType<float>
{
   /// \return 0
   constexpr static float zero () { return 0.0f; }
   /// \return 1
   constexpr static float one  () { return 1.0f; }

   template<class _Iterator>
   static void fill(_Iterator first, _Iterator last, const float& val)
   {
      static_assert(std::is_convertible<float, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      std::fill(first, last, val);
   }

   template<class _Iterator>
   static void scal(_Iterator first, _Iterator last, const float& val)
   {
      static_assert(std::is_convertible<float, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      while (first != last)
      {
         (*first) *= val;
         ++first;
      }
   }
};

/// Double precision real number
template <> struct NumericType<double>
{
   /// \return 0
   constexpr static double zero () { return 0.0; }
   /// \return 1
   constexpr static double one  () { return 1.0; }

   template<class _Iterator>
   static void fill(_Iterator first, _Iterator last, const double& val)
   {
      static_assert(std::is_convertible<double, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      std::fill(first, last, val);
   }

   template<class _Iterator>
   static void scal(_Iterator first, _Iterator last, const double& val)
   {
      static_assert(std::is_convertible<double, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      while (first != last)
      {
         (*first) *= val;
         ++first;
      }
   }
};

/// Single precision complex number
template <> struct NumericType<std::complex<float>>
{
   /// \return 0
   const static std::complex<float> zero () { return std::complex<float>(0.0, 0.0); }
   /// \return 1
   const static std::complex<float> one  () { return std::complex<float>(1.0, 0.0); }
   /// \return 1i
   const static std::complex<float> onei () { return std::complex<float>(0.0, 1.0); }

   template<class _Iterator>
   static void fill(_Iterator first, _Iterator last, const std::complex<float>& val)
   {
      static_assert(std::is_convertible<std::complex<float>, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      std::fill(first, last, val);
   }

   template<class _Iterator>
   static void scal(_Iterator first, _Iterator last, const std::complex<float>& val)
   {
      static_assert(std::is_convertible<std::complex<float>, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      while (first != last)
      {
         (*first) *= val;
         ++first;
      }
   }
};

/// Double precision complex number
template <> struct NumericType<std::complex<double>>
{
   /// \return 0
   const static std::complex<double> zero () { return std::complex<double>(0.0, 0.0); }
   /// \return 1
   const static std::complex<double> one  () { return std::complex<double>(1.0, 0.0); }
   /// \return 1i
   const static std::complex<double> onei () { return std::complex<double>(0.0, 1.0); }

   template<class _Iterator>
   static void fill(_Iterator first, _Iterator last, const std::complex<double>& val)
   {
      static_assert(std::is_convertible<std::complex<double>, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      std::fill(first, last, val);
   }

   template<class _Iterator>
   static void scal(_Iterator first, _Iterator last, const std::complex<double>& val)
   {
      static_assert(std::is_convertible<std::complex<double>, typename std::iterator_traits<_Iterator>::value_type>::value, "Value type is not convertible");
      while (first != last)
      {
         (*first) *= val;
         ++first;
      }
   }
};

}; // namespace btas

#endif // __BTAS_NUMERIC_TYPE_H
