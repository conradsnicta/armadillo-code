// Copyright (C) 2011-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_resize
//! @{



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_resize>
resize(const Base<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_resize>(X.get_ref(), in_n_rows, in_n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_resize>
resize(const Base<typename T1::elem_type,T1>& X, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_resize>(X.get_ref(), s.n_rows, s.n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const OpCube<T1, op_resize>
resize(const BaseCube<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_resize>(X.get_ref(), in_n_rows, in_n_cols, in_n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const OpCube<T1, op_resize>
resize(const BaseCube<typename T1::elem_type,T1>& X, const SizeCube& s)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_resize>(X.get_ref(), s.n_rows, s.n_cols, s.n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_resize>
resize(const SpBase<typename T1::elem_type, T1>& X, const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_resize>(X.get_ref(), in_n_rows, in_n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_resize>
resize(const SpBase<typename T1::elem_type, T1>& X, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_resize>(X.get_ref(), s.n_rows, s.n_cols);
  }



//! @}
