// Copyright (C) 2012-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup spop_misc
//! @{


class spop_scalar_times
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_scalar_times>& in);
  };



class spop_square
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_square>& in);
  };



class spop_sqrt
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sqrt>& in);
  };



class spop_abs
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_abs>& in);
  };



class spop_cx_abs
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>& in);
  };



class spop_arg
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_arg>& in);
  };



class spop_cx_arg
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_arg>& in);
  };



class spop_real
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_real>& in);
  };



class spop_imag
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_imag>& in);
  };



class spop_conj
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_conj>& in);
  };



class spop_repmat
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_repmat>& in);
  };



class spop_reshape
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_reshape>& in);
  };



class spop_resize
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_resize>& in);
  };



class spop_floor
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_floor>& in);
  };



class spop_ceil
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_ceil>& in);
  };



class spop_round
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_round>& in);
  };



class spop_trunc
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_trunc>& in);
  };



class spop_sign
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sign>& in);
  };



//! @}
