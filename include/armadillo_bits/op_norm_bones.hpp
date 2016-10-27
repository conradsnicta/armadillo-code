// Copyright (C) 2008-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_norm
//! @{


class op_norm
  {
  public:
  
  // norms for dense vectors and matrices
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_1(const Proxy<T1>& P, const typename  arma_not_cx<typename T1::elem_type>::result* junk = 0);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_1(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  template<typename eT> arma_hot inline static eT                    vec_norm_1_direct_std(const Mat<eT>& X);
  template<typename eT> arma_hot inline static eT                    vec_norm_1_direct_mem(const uword N, const eT* A);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_2(const Proxy<T1>& P, const typename  arma_not_cx<typename T1::elem_type>::result* junk = 0);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_2(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_std(const Mat<eT>& X);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_mem(const uword N, const eT* A);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_robust(const Mat<eT>& X);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_k(const Proxy<T1>& P, const int k);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_max(const Proxy<T1>& P);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_min(const Proxy<T1>& P);
  
  template<typename T1> inline static typename T1::pod_type mat_norm_1(const Proxy<T1>& P);
  template<typename T1> inline static typename T1::pod_type mat_norm_2(const Proxy<T1>& P);
  
  template<typename T1> inline static typename T1::pod_type mat_norm_inf(const Proxy<T1>& P);
  
  
  // norms for sparse matrices
  
  template<typename T1> inline static typename T1::pod_type mat_norm_1(const SpProxy<T1>& P);

  template<typename T1> inline static typename T1::pod_type mat_norm_2(const SpProxy<T1>& P, const typename arma_real_only<typename T1::elem_type>::result* junk = 0);
  template<typename T1> inline static typename T1::pod_type mat_norm_2(const SpProxy<T1>& P, const typename   arma_cx_only<typename T1::elem_type>::result* junk = 0);

  template<typename T1> inline static typename T1::pod_type mat_norm_inf(const SpProxy<T1>& P);
  };



//! @}
