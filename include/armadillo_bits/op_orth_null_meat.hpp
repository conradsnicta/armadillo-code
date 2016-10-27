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



template<typename T1>
inline
void
op_orth::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_orth>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const T tol = access::tmp_real(expr.aux);
  
  const bool status = op_orth::apply_direct(out, expr.m, tol);
  
  if(status == false)
    {
    arma_stop_runtime_error("orth(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_orth::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  arma_debug_check((tol < T(0)), "orth(): tolerance must be >= 0");
  
  const unwrap<T1>   tmp(expr.get_ref());
  const Mat<eT>& X = tmp.M;
  
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  const bool status = auxlib::svd_dc(U, s, V, X);
  
  V.reset();
  
  if(status == false)  { out.reset(); return false; }
  
  if(s.is_empty())  { out.reset(); return true; }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(X.n_rows, X.n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  if(count > 0)
    {
    out = U.head_cols(count);  // out *= eT(-1);
    }
  else
    {
    out.set_size(X.n_rows, 0);
    }
  
  return true;
  }



//



template<typename T1>
inline
void
op_null::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_null>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const T tol = access::tmp_real(expr.aux);
  
  const bool status = op_null::apply_direct(out, expr.m, tol);
  
  if(status == false)
    {
    arma_stop_runtime_error("null(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_null::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  arma_debug_check((tol < T(0)), "null(): tolerance must be >= 0");
  
  const unwrap<T1>   tmp(expr.get_ref());
  const Mat<eT>& X = tmp.M;
  
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  const bool status = auxlib::svd_dc(U, s, V, X);
  
  U.reset();
  
  if(status == false)  { out.reset(); return false; }
  
  if(s.is_empty())  { out.reset(); return true; }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(X.n_rows, X.n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  if(count < X.n_cols)
    {
    out = V.tail_cols(X.n_cols - count);
    
    const uword out_n_elem = out.n_elem;
          eT*   out_mem    = out.memptr();
    
    for(uword i=0; i<out_n_elem; ++i)
      {
      if(std::abs(out_mem[i]) < std::numeric_limits<T>::epsilon())  { out_mem[i] = eT(0); }
      }
    }
  else
    {
    out.set_size(X.n_cols, 0);
    }
  
  return true;
  }



//! @}
