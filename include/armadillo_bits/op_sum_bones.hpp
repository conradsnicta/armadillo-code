// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_sum
//! @{


class op_sum
  {
  public:
  
  // dense matrices
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_sum>& in);
  
  template<typename T1>
  arma_hot inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_unwrap(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  
  // cubes
  
  template<typename T1>
  arma_hot inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1, op_sum>& in);
  
  template<typename T1>
  arma_hot inline static void apply_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_unwrap(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_proxy(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const uword dim);
  };


//! @}
