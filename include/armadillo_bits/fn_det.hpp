// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_det
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, typename T1::elem_type >::result
det
  (
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return auxlib::det(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
det
  (
  const Op<T1, op_diagmat>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  arma_debug_check( (A.n_rows != A.n_cols), "det(): given matrix must be square sized" );
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  eT val1 = eT(1);
  eT val2 = eT(1);
  
  uword i,j;
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    val1 *= A[i];
    val2 *= A[j];
    }
  
  
  if(i < N)
    {
    val1 *= A[i];
    }
  
  return val1 * val2;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
det
  (
  const Op<T1, op_trimat>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.m);
  
  const uword N = P.get_n_rows();
  
  arma_debug_check( (N != P.get_n_cols()), "det(): given matrix must be square sized" );
  
  eT val1 = eT(1);
  eT val2 = eT(1);
  
  uword i,j;
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    val1 *= P.at(i,i);
    val2 *= P.at(j,j);
    }
  
  if(i < N)
    {
    val1 *= P.at(i,i);
    }
  
  return val1 * val2;
  }



//! determinant of inv(A), without doing the inverse operation
template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, typename T1::elem_type >::result
det
  (
  const Op<T1,op_inv>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT tmp = det(X.m);
  
  if(tmp == eT(0))  { arma_debug_warn("det(): denominator is zero" ); }
  
  return eT(1) / tmp;
  }



template<typename T1>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, typename T1::elem_type >::result
det
  (
  const Base<typename T1::elem_type,T1>& X,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("det(X,bool) is deprecated and will be removed; change to det(X)");
  
  return det(X.get_ref());
  }



template<typename T1>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, typename T1::elem_type >::result
det
  (
  const Base<typename T1::elem_type,T1>& X,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_warn("det(X,char*) is deprecated and will be removed; change to det(X)");
  
  return det(X.get_ref());
  }



template<typename T>
arma_warn_unused
arma_inline
const typename arma_scalar_only<T>::result &
det(const T& x)
  {
  return x;
  }



//! @}
