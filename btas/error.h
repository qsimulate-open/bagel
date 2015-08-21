#ifndef __BTAS_ERROR_H
#define __BTAS_ERROR_H

#include <exception>

namespace btas {

  /// exception class, used to mark exceptions specific to BTAS
  class exception : public std::exception {
    public:
      exception(const char* m) : message_(m) { }

      virtual const char* what() const noexcept { return message_; }

    private:
      const char* message_;
  }; // class exception

  /// Place a break point in this function to stop before btas::exception is thrown.
  inline void exception_break() { }

} // namespace btas

// configure BTAS_ASSERT
#ifdef BTAS_ASSERT_THROWS

#define BTAS_STRINGIZE( s ) #s

#define BTAS_EXCEPTION_MESSAGE( file , line , mess ) \
  "BTAS: exception at " file "(" BTAS_STRINGIZE( line ) "): " mess ". Break in btas::exception_break to learn more."

#define BTAS_EXCEPTION( m ) \
    { \
      btas::exception_break(); \
      throw btas::exception ( BTAS_EXCEPTION_MESSAGE( __FILE__ , __LINE__ , m ) ); \
    }

#define BTAS_ASSERT( a )  if(! ( a ) ) BTAS_EXCEPTION( "assertion failed" )

#else // defined BTAS_ASSERT_THROWS

#define BTAS_ASSERT( a )  assert((a));

#endif // defined BTAS_ASSERT_THROWS

#endif // __BTAS_ERROR_H
