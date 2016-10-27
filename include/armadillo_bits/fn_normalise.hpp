// Copyright (C) 2014-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_normalise
//! @{



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (resolves_to_vector<T1>::value == true)),
  const Op<T1, op_normalise_vec>
  >::result
normalise
  (
  const T1&   X,
  const uword p = uword(2),
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return Op<T1, op_normalise_vec>(X, p, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (resolves_to_vector<T1>::value == false)),
  const Op<T1, op_normalise_mat>
  >::result
normalise
  (
  const T1&   X,
  const uword p = uword(2),
  const uword dim = 0,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_normalise_mat>(X, p, dim);
  }



//! @}
