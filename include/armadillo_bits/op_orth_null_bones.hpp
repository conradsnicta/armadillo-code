// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup op_orth_null
//! @{



class op_orth
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_orth>& expr);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol);
  };



class op_null
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_null>& expr);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol);
  };



//! @}
