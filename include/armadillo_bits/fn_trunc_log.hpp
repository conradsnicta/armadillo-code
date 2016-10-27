// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_trunc_log
//! @{



template<typename eT>
arma_warn_unused
inline
static
typename arma_real_only<eT>::result
trunc_log(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559)
    {
    if(x == std::numeric_limits<eT>::infinity())
      {
      return Datum<eT>::log_max;
      }
    else
      {
      return (x <= eT(0)) ? Datum<eT>::log_min : std::log(x);
      }
    }
  else
    {
    return std::log(x);
    }
  }



template<typename eT>
arma_warn_unused
inline
static
typename arma_integral_only<eT>::result
trunc_log(const eT x)
  {
  return eT( trunc_log( double(x) ) );
  }



template<typename T>
arma_warn_unused
inline
static
std::complex<T>
trunc_log(const std::complex<T>& x)
  {
  return std::complex<T>( trunc_log( std::abs(x) ), std::arg(x) );
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOp<T1, eop_trunc_log>
trunc_log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_log>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_trunc_log>
trunc_log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_trunc_log>(A.get_ref());
  }



//! @}
