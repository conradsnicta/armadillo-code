// Copyright (C) 2009-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Dimitrios Bouzas



//! \addtogroup fn_shuffle
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value),
  const Op<T1, op_shuffle_default>
  >::result
shuffle
  (
  const T1& X
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_shuffle_default>(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value),
  const Op<T1, op_shuffle>
  >::result
shuffle
  (
  const T1&   X,
  const uword dim
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_shuffle>(X, dim, 0);
  }



//! @}
