// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_logmat
//! @{



class op_logmat
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat< std::complex<typename T1::elem_type> >& out, const mtOp<std::complex<typename T1::elem_type>,T1,op_logmat>& in);

  template<typename T1>
  inline static bool apply_direct(Mat< std::complex<typename T1::elem_type> >& out, const Op<T1,op_diagmat>& expr, const uword);
  
  template<typename T1>
  inline static bool apply_direct(Mat< std::complex<typename T1::elem_type> >& out, const Base<typename T1::elem_type,T1>& expr, const uword n_iters);
  };



class op_logmat_cx
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_logmat_cx>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& expr, const uword);
  
  template<typename T1>
  inline static bool apply_direct_noalias(Mat<typename T1::elem_type>& out, const diagmat_proxy<T1>& P);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const uword n_iters);
  
  template<typename T>
  inline static bool apply_common(Mat< std::complex<T> >& out, Mat< std::complex<T> >& S, const uword n_iters);
  
  
  template<typename eT>
  inline static bool helper(Mat<eT>& S, const uword m);
  };



class op_logmat_sympd
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_logmat_sympd>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr);
  };



//! @}
