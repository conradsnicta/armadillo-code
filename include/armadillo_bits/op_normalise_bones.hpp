// Copyright (C) 2014-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup op_normalise
//! @{



class op_normalise_vec
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_normalise_vec>& in);
  };



class op_normalise_mat
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_normalise_mat>& in);
  
  template<typename eT> inline static void apply(Mat<eT>& out, const Mat<eT>& A, const uword p, const uword dim);
  };



//! @}
