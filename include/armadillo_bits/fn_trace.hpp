// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup fn_trace
//! @{


template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.get_ref());
  
  const uword N = (std::min)(P.get_n_rows(), P.get_n_cols());
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  uword i,j;
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    val1 += P.at(i,i);
    val2 += P.at(j,j);
    }
  
  if(i < N)
    {
    val1 += P.at(i,i);
    }
  
  return val1 + val2;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  eT val = eT(0);
  
  for(uword i=0; i<N; ++i)
    {
    val += A[i];
    }
  
  return val;
  }



//! speedup for trace(A*B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> tmp1(X.A);
  const partial_unwrap<T2> tmp2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = tmp1.M;
  const typename partial_unwrap<T2>::stored_type& B = tmp2.M;
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )
    {
    return eT(0);
    }
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;

  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  eT acc = eT(0);
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_cols = B_n_rows
      
      uword j;
      
      for(j=1; j < A_n_cols; j+=2)
        {
        const uword i = (j-1);
        
        const eT tmp_i = B_colptr[i];
        const eT tmp_j = B_colptr[j];
        
        acc1 += A.at(k, i) * tmp_i;
        acc2 += A.at(k, j) * tmp_j;
        }
      
      const uword i = (j-1);
      
      if(i < A_n_cols)
        {
        acc1 += A.at(k, i) * B_colptr[i];
        }
      }
      
    acc = (acc1 + acc2);
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_rows = B_n_rows
      acc += op_dot::direct_dot(A_n_rows, A_colptr, B_colptr);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      // condition: A_n_cols = B_n_cols
      for(uword i=0; i < A_n_cols; ++i)
        {
        acc += A.at(k,i) * B.at(k,i);
        }
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      
      // condition: A_n_rows = B_n_cols
      for(uword i=0; i < A_n_rows; ++i)
        {
        acc += A_colptr[i] * B.at(k,i);
        }
      }
    }
  
  return (use_alpha) ? (alpha * acc) : acc;
  }



//! trace of sparse object; generic version
template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpBase<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> P(expr.get_ref());
  
  const uword N = (std::min)(P.get_n_rows(), P.get_n_cols());
  
  eT acc = eT(0);
  
  if( (is_SpMat<typename SpProxy<T1>::stored_type>::value) && (P.get_n_nonzero() >= 5*N) )
    {
    const unwrap_spmat<typename SpProxy<T1>::stored_type> U(P.Q);
    
    const SpMat<eT>& X = U.M;
    
    for(uword i=0; i < N; ++i)
      {
      acc += X.at(i,i);  // use binary search
      }
    }
  else
    {
    typename SpProxy<T1>::const_iterator_type it = P.begin();
    
    const uword P_n_nz = P.get_n_nonzero();
    
    for(uword i=0; i < P_n_nz; ++i)
      {
      if(it.row() == it.col())  { acc += (*it); }
      
      ++it;
      }
    }
  
  return acc;
  }



//! trace of sparse object; speedup for trace(A*B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpGlue<T1, T2, spglue_times>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // better-than-nothing implementation
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_nonzero == 0) || (B.n_nonzero == 0) )
    {
    return eT(0);
    }
  
  const uword N = (std::min)(A.n_rows, B.n_cols);
  
  eT acc = eT(0);
  
  for(uword k=0; k < N; ++k)
    {
    typename SpMat<eT>::const_col_iterator B_it     = B.begin_col(k);
    typename SpMat<eT>::const_col_iterator B_it_end = B.end_col(k);
    
    while(B_it != B_it_end)
      {
      const eT    B_val = (*B_it);
      const uword i     = B_it.row();
      
      acc += A.at(k,i) * B_val;
      
      ++B_it;
      }
    }
  
  return acc;
  }



//! @}
