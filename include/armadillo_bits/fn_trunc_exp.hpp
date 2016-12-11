// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_trunc_exp
//! @{



template<typename eT>
arma_warn_unused
inline
static
typename arma_real_only<eT>::result
trunc_exp(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559 && (x >= Datum<eT>::log_max ))
    {
    return std::numeric_limits<eT>::max();
    }
  else
    {
    return std::exp(x);
    }
  }



template<typename eT>
arma_warn_unused
inline
static
typename arma_integral_only<eT>::result
trunc_exp(const eT x)
  {
  return eT( trunc_exp( double(x) ) );
  }



template<typename T>
arma_warn_unused
inline
static
std::complex<T>
trunc_exp(const std::complex<T>& x)
  {
  return std::polar( trunc_exp( x.real() ), x.imag() );
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_trunc_exp> >::result
trunc_exp(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_exp>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_trunc_exp>
trunc_exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_trunc_exp>(A.get_ref());
  }



//! @}
