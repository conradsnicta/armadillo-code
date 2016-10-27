// Copyright (C) 2008-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_inv
//! @{


//! immediate inverse of a matrix, storing the result in a dense matrix
template<typename eT>
inline
void
op_inv::apply(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  // no need to check for aliasing, due to:
  // - auxlib::inv() copies A to out before inversion
  // - for 2x2 and 3x3 matrices the code is alias safe
  
  bool status = auxlib::inv(out, A);
  
  if(status == false)
    {
    out.reset();
    arma_stop_runtime_error("inv(): matrix seems singular");
    }
  }



//! immediate inverse of T1, storing the result in a dense matrix
template<typename T1>
inline
void
op_inv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& X)
  {
  arma_extra_debug_sigprint();
  
  const strip_diagmat<T1> strip(X.m);
  
  bool status;
  
  if(strip.do_diagmat == true)
    {
    status = op_inv::apply_diagmat(out, strip.M);
    }
  else
    {
    status = auxlib::inv(out, X.m);
    }
    
  if(status == false)
    {
    out.reset();
    arma_stop_runtime_error("inv(): matrix seems singular");
    }
  }



template<typename T1>
inline
bool
op_inv::apply_diagmat(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X);
  
  arma_debug_check( (A.n_rows != A.n_cols), "inv(): given matrix must be square sized" );
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  bool status = true;
  
  if(A.is_alias(out) == false)
    {
    out.zeros(N,N);
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = A[i];
      
      out.at(i,i) = eT(1) / val;
      
      if(val == eT(0))  { status = false; }
      }
    }
  else
    {
    Mat<eT> tmp(N, N, fill::zeros);
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = A[i];
      
      tmp.at(i,i) = eT(1) / val;
      
      if(val == eT(0))  { status = false; }
      }
    
    out.steal_mem(tmp);
    }
  
  return status;
  }



//! inverse of T1 (triangular matrices)
template<typename T1>
inline
void
op_inv_tr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_tr>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = auxlib::inv_tr(out, X.m, X.aux_uword_a);
  
  if(status == false)
    {
    out.reset();
    arma_stop_runtime_error("inv(): matrix seems singular");
    }
  }



//! inverse of T1 (symmetric positive definite matrices)
template<typename T1>
inline
void
op_inv_sympd::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_sympd>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = auxlib::inv_sympd(out, X.m);
  
  if(status == false)
    {
    out.reset();
    arma_stop_runtime_error("inv_sympd(): matrix is singular or not positive definite");
    }
  }



//! @}
