// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_find_unique
//! @{



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<uword,T1,op_find_unique>
  >::result
find_unique
  (
  const T1&  X,
  const bool ascending_indices = true
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword,T1,op_find_unique>(X, ((ascending_indices) ? uword(1) : uword(0)), uword(0));
  }



template<typename T1>
arma_warn_unused
inline
uvec
find_unique
  (
  const BaseCube<typename T1::elem_type,T1>& X,
  const bool ascending_indices = true
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X.get_ref());
  
  const Mat<eT> R( const_cast< eT* >(tmp.M.memptr()), tmp.M.n_elem, 1, false );
  
  return find_unique(R,ascending_indices);
  }



//! @}
