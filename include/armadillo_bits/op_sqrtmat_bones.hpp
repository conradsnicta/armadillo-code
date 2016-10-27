// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_sqrtmat
//! @{



class op_sqrtmat
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat< std::complex<typename T1::elem_type> >& out, const mtOp<std::complex<typename T1::elem_type>,T1,op_sqrtmat>& in);

  template<typename T1>
  inline static bool apply_direct(Mat< std::complex<typename T1::elem_type> >& out, const Op<T1,op_diagmat>& expr);
  
  template<typename T1>
  inline static bool apply_direct(Mat< std::complex<typename T1::elem_type> >& out, const Base<typename T1::elem_type,T1>& expr);
  };



class op_sqrtmat_cx
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sqrtmat_cx>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& expr);
  
  template<typename T1>
  inline static bool apply_direct_noalias(Mat<typename T1::elem_type>& out, const diagmat_proxy<T1>& P);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr);
  
  template<typename T>
  inline static bool helper(Mat< std::complex<T> >& S);
  };



class op_sqrtmat_sympd
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sqrtmat_sympd>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr);
  };



//! @}
