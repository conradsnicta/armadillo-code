// Copyright (C) 2010-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_conv
//! @{



//! Convolution, which is also equivalent to polynomial multiplication and FIR digital filtering.

template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, T2, glue_conv>
  >::result
conv(const T1& A, const T2& B, const char* shape = "full")
  {
  arma_extra_debug_sigprint();
  
  const char sig = (shape != NULL) ? shape[0] : char(0);
  
  arma_debug_check( ((sig != 'f') && (sig != 's')), "conv(): unsupported value of 'shape' parameter" );
  
  const uword mode = (sig == 's') ? uword(1) : uword(0);
  
  return Glue<T1, T2, glue_conv>(A, B, mode);
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, T2, glue_conv2>
  >::result
conv2(const T1& A, const T2& B, const char* shape = "full")
  {
  arma_extra_debug_sigprint();
  
  const char sig = (shape != NULL) ? shape[0] : char(0);
  
  arma_debug_check( ((sig != 'f') && (sig != 's')), "conv2(): unsupported value of 'shape' parameter" );
  
  const uword mode = (sig == 's') ? uword(1) : uword(0);
  
  return Glue<T1, T2, glue_conv2>(A, B, mode);
  }



//! @}
