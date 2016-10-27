// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_diff
//! @{



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_diff_default>
  >::result
diff
  (
  const T1&   X,
  const uword k = 1
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diff_default>(X, k, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_diff>
  >::result
diff
  (
  const T1&   X,
  const uword k,
  const uword dim
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diff>(X, k, dim);
  }



//! @}
