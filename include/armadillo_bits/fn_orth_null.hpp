// Copyright (C) 2015-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_orth_null
//! @{



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_real<typename T1::pod_type>::value, const Op<T1, op_orth> >::result
orth(const Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tol = 0.0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return Op<T1, op_orth>(X.get_ref(), eT(tol));
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
orth(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tol = 0.0)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_orth::apply_direct(out, X.get_ref(), tol);
  
  if(status == false)
    {
    arma_debug_warn("orth(): svd failed");
    }
  
  return status;
  }



//



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_real<typename T1::pod_type>::value, const Op<T1, op_null> >::result
null(const Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tol = 0.0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return Op<T1, op_null>(X.get_ref(), eT(tol));
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
null(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tol = 0.0)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_null::apply_direct(out, X.get_ref(), tol);
  
  if(status == false)
    {
    arma_debug_warn("null(): svd failed");
    }
  
  return status;
  }



//! @}
