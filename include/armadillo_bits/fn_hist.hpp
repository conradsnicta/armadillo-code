// Copyright (C) 2012-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_hist
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value) && (is_not_complex<typename T1::elem_type>::value),
  const mtOp<uword,T1,op_hist>
  >::result
hist(const T1& A, const uword n_bins = 10)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword,T1,op_hist>(A, n_bins, 0);
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value) && (is_arma_type<T2>::value) && (is_not_complex<typename T1::elem_type>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const mtGlue<uword,T1,T2,glue_hist_default>
  >::result
hist(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword,T1,T2,glue_hist_default>(X, Y);
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value) && (is_arma_type<T2>::value) && (is_not_complex<typename T1::elem_type>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const mtGlue<uword,T1,T2,glue_hist>
  >::result
hist(const T1& X, const T2& Y, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword,T1,T2,glue_hist>(X, Y, dim);
  }


//! @}
