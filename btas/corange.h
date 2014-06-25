/*
 * corange.h
 *
 *  Created on: Dec 27, 2013
 *      Author: evaleev
 */

#ifndef BTAS_CORANGE_H_
#define BTAS_CORANGE_H_

#include <tuple>
#include <btas/util/functional.h>

namespace btas {

  /// CoRange is a pack of Range objects
  template <typename ...Ranges>
  class CoRange;

  /// CoRangeIterator implements iteration over CoRange
  /// it resembles Boost.Iterator's <a href="http://www.boost.org/doc/libs/1_55_0/libs/iterator/doc/zip_iterator.html">zip_iterator</a>.
  template <typename ...Ranges>
  class CoRangeIterator;

  /// CoRangeIterator over a pair of ranges.
  template <typename R1, typename R2>
  class CoRangeIterator<R1,R2> {
    public:
      typedef CoRange<R1, R2> corange_type;
      typedef std::tuple<typename R1::iterator, typename R2::iterator> iterator;
      typedef std::tuple< typename R1::iterator::value_type, typename R2::iterator::value_type> value_type;

      CoRangeIterator(corange_type& corange,
                      const iterator& iter) : corange_(corange), iter_(iter) {}

      value_type operator*() const {
        return std::make_tuple(*(std::get<0>(iter_)), *(std::get<1>(iter_)));
      }

      iterator iter() {
        return iter_;
      }

      /// reset second iterator to begin, if reached end
      void operator++() {
        ++(std::get<0>(iter_));
        ++(std::get<1>(iter_));
        if (std::get<1>(iter_) == std::get<1>(corange_.ranges()).end())
          std::get<1>(iter_) = std::get<1>(corange_.ranges()).begin();
      }

      template <typename _R1, typename _R2>
      friend bool operator==(const CoRangeIterator<_R1,_R2>& i1,
          const CoRangeIterator<_R1,_R2>& i2);

    private:
      corange_type& corange_;
      iterator iter_;
  };

  template <typename R1, typename R2>
  bool operator==(const CoRangeIterator<R1,R2>& i1,
                  const CoRangeIterator<R1,R2>& i2) {
    return std::get<0>(i1.iter_) == std::get<0>(i2.iter_);
  }
  template <typename R1, typename R2>
  bool operator!=(const CoRangeIterator<R1,R2>& i1,
                  const CoRangeIterator<R1,R2>& i2) {
    return not operator==(i1,i2);
  }

  /// This is a CoRange of two Ranges. The first Range iterates once from begin to end; the second Range
  /// re-starts from begin if necessary.
  template <typename R1, typename R2>
  class CoRange<R1,R2> {
    public:

      typedef std::tuple<R1&,R2&> ranges_type;

      typedef CoRangeIterator<R1,R2> iterator;
      typedef CoRangeIterator<const R1,const R2> const_iterator;

      CoRange(R1& r1, R2& r2) : ranges_(r1,r2) {}

      const ranges_type& ranges() const { return ranges_; }

      iterator begin() {
        return iterator(*this, std::make_tuple(std::get<0>(ranges_).begin(), std::get<1>(ranges_).begin()));
      }

      const_iterator begin() const {
        return const_iterator(*this, std::make_tuple(std::get<0>(ranges_).begin(), std::get<1>(ranges_).begin()));
      }

      iterator end() {
        return iterator(*this, std::make_tuple(std::get<0>(ranges_).end(), std::get<1>(ranges_).end()));
      }

      const_iterator end() const {
        return const_iterator(*this, std::make_tuple(std::get<0>(ranges_).end(), std::get<1>(ranges_).end()));
      }

    private:
      ranges_type ranges_;
  };

  template <typename ... Ranges>
  CoRange<Ranges...> make_corange(Ranges&... args) {
    return CoRange<Ranges...>(args...);
  }


#if 0   /// This is hard or impossible ... need to make a tuple by applying a function to a tuple

  /// CoRangeIterator iterates over the pack of Ranges
  template <typename ...Ranges>
    class CoRangeIterator : public std::tuple<typename Ranges::iterator ...> {
    public:
      typedef CoRange<Ranges...> corange_type;

      typedef std::tuple< typename Ranges::iterator::value_type ...> value_type;

      // how???
      value_type operator*() const {
        return
      }
  };

  /// CoRange is a pack of Range objects
  template <typename ...Ranges>
  class CoRange {
    public:

      typedef std::tuple<Ranges&...> ranges_type;

      typedef CoRangeIterator<Ranges...> iterator;

      CoRange(Ranges&... ranges) : ranges_(ranges...) {}

      void print(std::ostream& os) {
        for(auto i=0; i<std::tuple_size<ranges_type>::value;
            ++i) {
          os << "CoRange[" << i << "]: ";
          switch (i) {
            case 0:
              os << std::get<0>(ranges_) << std::endl; break;
            case 1:
              os << std::get<1>(ranges_) << std::endl; break;
            default:
              assert(false);
          }
        }
      }

    private:
      std::tuple<Ranges&...> ranges_;
  };
#endif

}


#endif /* BTAS_CORANGE_H_ */
