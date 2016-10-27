// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_diagmat
//! @{


//! interpret a matrix or a vector as a diagonal matrix (i.e. off-diagonal entries are zero)
template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_diagmat>
  >::result
diagmat(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diagmat>(X);
  }



//! create a matrix with the k-th diagonal set to the given vector
template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_diagmat2>
  >::result
diagmat(const T1& X, const sword k)
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = (k < 0) ? uword(-k) : uword(0);
  const uword col_offset = (k > 0) ? uword( k) : uword(0);
  
  return Op<T1, op_diagmat2>(X, row_offset, col_offset);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_diagmat>
diagmat(const SpBase<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_diagmat>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_diagmat2>
diagmat(const SpBase<typename T1::elem_type,T1>& X, const sword k)
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = (k < 0) ? uword(-k) : uword(0);
  const uword col_offset = (k > 0) ? uword( k) : uword(0);
  
  return SpOp<T1, spop_diagmat2>(X.get_ref(), row_offset, col_offset);
  }



//! @}
