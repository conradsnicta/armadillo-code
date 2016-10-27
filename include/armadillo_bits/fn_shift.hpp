// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup fn_shift
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value),
  const Op<T1, op_shift_default>
  >::result
shift
  (
  const T1&   X,
  const sword N
  )
  {
  arma_extra_debug_sigprint();
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  return Op<T1, op_shift_default>(X, len, neg, uword(0), 'j');
  }



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value),
  const Op<T1, op_shift>
  >::result
shift
  (
  const T1&   X,
  const sword N,
  const uword dim
  )
  {
  arma_extra_debug_sigprint();
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  return Op<T1, op_shift>(X, len, neg, dim, 'j');
  }



//! @}
