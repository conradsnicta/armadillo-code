// Copyright (C) 2010-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup glue_join
//! @{



class glue_join_cols
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_join_cols>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& A, const Proxy<T2>& B);
  };



class glue_join_rows
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_join_rows>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& A, const Proxy<T2>& B);
  };



class glue_join_slices
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_join_slices>& X);
  };



//! @}

