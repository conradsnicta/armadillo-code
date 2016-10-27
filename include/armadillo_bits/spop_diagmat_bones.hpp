// Copyright (C) 2012-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup spop_diagmat
//! @{


class spop_diagmat
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_diagmat>& in);
  
  template<typename T1>
  inline static void apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p);
  };



class spop_diagmat2
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_diagmat2>& in);
  
  template<typename eT>
  inline static void apply_noalias(SpMat<eT>& out, const SpMat<eT>& X, const uword row_offset, const uword col_offset);
  };



//! @}
