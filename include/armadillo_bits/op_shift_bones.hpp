// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup op_shift
//! @{



class op_shift_default
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_shift_default>& in);
  };



class op_shift
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_shift>& in);
  
  template<typename eT> inline static void apply_direct(Mat<eT>& out, const Mat<eT>& X, const uword len, const uword neg, const uword dim);
  
  template<typename eT> inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword len, const uword neg, const uword dim);
  
  template<typename eT> inline static void apply_alias(Mat<eT>& out, const uword len, const uword neg, const uword dim);
  };



//! @}
