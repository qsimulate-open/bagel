#ifndef __BTAS_BTAS_H
#define __BTAS_BTAS_H

#include <cassert>

#include <btas/generic/dot_impl.h>
#include <btas/generic/scal_impl.h>
#include <btas/generic/axpy_impl.h>
#include <btas/generic/ger_impl.h>
#include <btas/generic/gemv_impl.h>
#include <btas/generic/gemm_impl.h>

#ifdef _CONTRACT_OPT_BAGEL
#include <btas/optimize/contract.h>
#else
#include <btas/generic/contract.h>
#endif

#endif // __BTAS_BTAS_H
