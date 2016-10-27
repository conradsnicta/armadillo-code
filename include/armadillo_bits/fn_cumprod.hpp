// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_cumprod
//! @{



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_cumprod_default>
  >::result
cumprod(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_cumprod_default>(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_cumprod>
  >::result
cumprod(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_cumprod>(X, dim, 0);
  }



template<typename T>
arma_warn_unused
arma_inline
const typename arma_scalar_only<T>::result &
cumprod(const T& x)
  {
  return x;
  }



//! @}
