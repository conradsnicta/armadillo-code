// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_index_min
//! @{


class op_index_min
  {
  public:
  
  // dense matrices
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword,T1,op_index_min>& in);
  
  template<typename eT>
  inline static void apply_noalias(Mat<uword>& out, const Mat<eT>& X, const uword dim);
  
  
  // sparse matrices
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const SpBase<typename T1::elem_type,T1>& expr, const uword dim);
  };



//! @}
