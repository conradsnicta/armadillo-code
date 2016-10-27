// Copyright (C) 2010-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup glue_conv
//! @{



class glue_conv
  {
  public:

  template<typename eT> inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const bool A_is_col);
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv>& X);
  };



class glue_conv2
  {
  public:
  
  template<typename eT> inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv2>& expr);
  };



//! @}

