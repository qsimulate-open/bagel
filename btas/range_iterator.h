/*
 * This file is a part of TiledArray.
 * Copyright (C) 2013  Virginia Tech
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef BTAS_RANGE_ITERATOR_H__INCLUDED
#define BTAS_RANGE_ITERATOR_H__INCLUDED

#include <iterator>

namespace btas {

    template <typename, typename>
    class RangeIterator;

} // namespace btas

namespace std {
  template <typename Value, typename Range>
  void advance(btas::RangeIterator<Value, Range>&,
      typename btas::RangeIterator<Value, Range>::difference_type );

  template <typename Value, typename Range>
  typename btas::RangeIterator<Value, Range>::difference_type
  distance(const btas::RangeIterator<Value, Range>&,
      const btas::RangeIterator<Value, Range>&);

} // namespace std

namespace btas {

    /// Iterates over a Range of Values

    /// This is an input iterator that is used to iterate over elements of a \c Range.
    /// \tparam Value The value type of the iterator
    /// \tparam Range The range that the iterator references
    /// \note The range object must define the function
    /// \c Range::increment(Value&) \c const, and be accessible to
    /// \c RangeIterator.
    template <typename Value,
              typename Range>
    class RangeIterator {
    public:

      // Standard iterator typedefs
      typedef Value value_type;                                    ///< Iterator value type
      typedef const value_type& reference;                         ///< Iterator reference type
      typedef const value_type* pointer;                           ///< Iterator pointer type
      typedef std::input_iterator_tag iterator_category;           ///< Iterator category tag
      typedef std::ptrdiff_t difference_type;                      ///< Iterator difference type

      // additional typedefs
      typedef Range range_type;                                    ///< Range type

      /// Copy constructor

      /// \param other The other iterator to be copied
      RangeIterator(const RangeIterator& other) :
        range_(other.range_), current_(other.current_) {
      }

      /// Construct an index iterator

      /// \param v The initial value of the iterator index
      /// \param c The range that the iterator will reference
      RangeIterator(const value_type& v, const range_type* c) :
          range_(c), current_(v)
      { }

      /// Copy constructor

      /// \param other The other iterator to be copied
      /// \return A reference to this object
      RangeIterator& operator=(const RangeIterator& other) {
        current_ = other.current_;
        range_ = other.range_;

        return *this;
      }

      const range_type* range() const { return range_; }

      /// Dereference operator

      /// \return A \c reference to the current data
      reference operator*() const { return current_; }

      /// Increment operator

      /// Increment the iterator
      /// \return The modified iterator
      RangeIterator& operator++() {
        range_->increment(current_);
        return *this;
      }

      /// Increment operator

      /// Increment the iterator
      /// \return An unmodified copy of the iterator
      RangeIterator operator++(int) {
        RangeIterator temp(*this);
        range_->increment(current_);
        return temp;
      }

      /// Pointer operator

      /// \return A \c pointer to the current data
      pointer operator->() const { return & current_; }

      void advance(difference_type n) {
        range_->advance(current_, n);
      }

      difference_type distance_to(const RangeIterator& other) const {
        TA_ASSERT(range_ == other.range_);
        return range_->distance_to(current_, other.current_);
      }

    private:

      const range_type* range_;          ///< The range that the iterator references
      value_type current_;               ///< The current value of the iterator
    }; // class RangeIterator

    /// Equality operator

    /// Compares the iterators for equality. They must reference the same range
    /// object to be considered equal.
    /// \tparam Value The value type of the iterator
    /// \tparam Range The range that the iterator references
    /// \param left_it The left-hand iterator to be compared
    /// \param right_it The right-hand iterator to be compared
    /// \return \c true if the the value and range are equal for the \c left_it
    /// and \c right_it , otherwise \c false .
    template <typename Value, typename Range>
    bool operator==(const RangeIterator<Value, Range>& left_it, const RangeIterator<Value, Range>& right_it) {
      return ((*left_it) == (*right_it)) &&
          (left_it.range() == right_it.range());
    }

    /// Inequality operator

    /// Compares the iterators for inequality.
    /// \tparam Value The value type of the iterator
    /// \tparam Range The range that the iterator references
    /// \param left_it The left-hand iterator to be compared
    /// \param right_it The right-hand iterator to be compared
    /// \return \c true if the the value or range are not equal for the
    /// \c left_it and \c right_it , otherwise \c false .
    template <typename Value, typename Range>
    bool operator!=(const RangeIterator<Value, Range>& left_it, const RangeIterator<Value, Range>& right_it) {
      return ((*left_it) != (*right_it)) ||
          (left_it.range() != right_it.range());
    }

} // namespace btas

namespace std {
  template <typename Value, typename Range>
  void advance(btas::RangeIterator<Value, Range>& it,
      typename btas::RangeIterator<Value, Range>::difference_type n)
  {
    it.advance(n);
  }

  template <typename Value, typename Range>
  typename btas::RangeIterator<Value, Range>::difference_type
  distance(const btas::RangeIterator<Value, Range>& first,
      const btas::RangeIterator<Value, Range>& last)
  {
    return first.distance_to(last);
  }

} // namespace std

#endif // BTAS_RANGE_ITERATOR_H__INCLUDED
