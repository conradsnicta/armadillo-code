// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_histc
//! @{


template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value) && (is_arma_type<T2>::value) && (is_not_complex<typename T1::elem_type>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const mtGlue<uword,T1,T2,glue_histc_default>
  >::result
histc(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword,T1,T2,glue_histc_default>(X, Y);
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value) && (is_arma_type<T2>::value) && (is_not_complex<typename T1::elem_type>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const mtGlue<uword,T1,T2,glue_histc>
  >::result
histc(const T1& X, const T2& Y, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword,T1,T2,glue_histc>(X, Y, dim);
  }


//! @}
