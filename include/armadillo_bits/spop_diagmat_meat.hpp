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


//! \addtogroup spop_diagmat
//! @{



template<typename T1>
inline
void
spop_diagmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    spop_diagmat::apply_noalias(out, p);
    }
  else
    {
    SpMat<eT> tmp;
    
    spop_diagmat::apply_noalias(tmp, p);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
inline
void
spop_diagmat::apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  const uword P_n_nz   = P.get_n_nonzero();
  
  const bool P_is_vec = (P_n_rows == 1) || (P_n_cols == 1);
  
  if(P_is_vec)    // generate a diagonal matrix out of a vector
    {
    const uword N = (P_n_rows == 1) ? P_n_cols : P_n_rows;
    
    out.zeros(N, N);
    
    if(P_n_nz == 0)  { return; }
    
    typename SpProxy<T1>::const_iterator_type it = P.begin();
    
    if(P_n_cols == 1)
      {
      for(uword i=0; i < P_n_nz; ++i)
        {
        const uword row = it.row();
        
        out.at(row,row) = (*it);
        
        ++it;
        }
      }
    else
    if(P_n_rows == 1)
      {
      for(uword i=0; i < P_n_nz; ++i)
        {
        const uword col = it.col();
        
        out.at(col,col) = (*it);
        
        ++it;
        }
      }
    }
  else   // generate a diagonal matrix out of a matrix
    {
    out.zeros(P_n_rows, P_n_cols);
    
    const uword N = (std::min)(P_n_rows, P_n_cols);
    
    if( (is_SpMat<typename SpProxy<T1>::stored_type>::value) && (P_n_nz >= 5*N) )
      {
      const unwrap_spmat<typename SpProxy<T1>::stored_type> U(P.Q);
      
      const SpMat<eT>& X = U.M;
      
      for(uword i=0; i < N; ++i)
        {
        const eT val = X.at(i,i);  // use binary search
        
        if(val != eT(0))  { out.at(i,i) = val; }
        }
      }
    else
      {
      if(P_n_nz == 0)  { return; }
      
      typename SpProxy<T1>::const_iterator_type it = P.begin();
      
      for(uword i=0; i < P_n_nz; ++i)
        {
        const uword row = it.row();
        const uword col = it.col();
        
        if(row == col)  { out.at(row,row) = (*it); }
        
        ++it;
        }
      }
    }
  }



template<typename T1>
inline
void
spop_diagmat2::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_diagmat2>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword row_offset = in.aux_uword_a;
  const uword col_offset = in.aux_uword_b;
  
  const unwrap_spmat<T1> U(in.m);
  
  if(&(U.M) == &out)
    {
    SpMat<eT> tmp;
    
    spop_diagmat2::apply_noalias(tmp, U.M, row_offset, col_offset);
    
    out.steal_mem(tmp);
    }
  else
    {
    spop_diagmat2::apply_noalias(out, U.M, row_offset, col_offset);
    }
  }



template<typename eT>
inline
void
spop_diagmat2::apply_noalias(SpMat<eT>& out, const SpMat<eT>& X, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = X.n_rows;
  const uword n_cols = X.n_cols;
  const uword n_elem = X.n_elem;
  
  if(n_elem == 0)  { out.reset(); return; }
  
  const bool X_is_vec = (n_rows == 1) || (n_cols == 1);
  
  if(X_is_vec)    // generate a diagonal matrix out of a vector
    {
    const uword n_pad = (std::max)(row_offset, col_offset);
    
    out.zeros(n_elem + n_pad, n_elem + n_pad);
    
    const uword X_n_nz = X.n_nonzero;
    
    if(X_n_nz == 0)  { return; }
    
    typename SpMat<eT>::const_iterator it = X.begin();
      
    if(n_cols == 1)
      {
      for(uword i=0; i < X_n_nz; ++i)
        {
        const uword row = it.row();
        
        out.at(row_offset + row, col_offset + row) = (*it);
        
        ++it;
        }
      }
    else
    if(n_rows == 1)
      {
      for(uword i=0; i < X_n_nz; ++i)
        {
        const uword col = it.col();
        
        out.at(row_offset + col, col_offset + col) = (*it);
        
        ++it;
        }
      }
    }
  else   // generate a diagonal matrix out of a matrix
    {
    arma_debug_check
      (
      ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
      "diagmat(): requested diagonal out of bounds"
      );
    
    out.zeros(n_rows, n_cols);
    
    if(X.n_nonzero == 0)  { return; }
    
    const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    for(uword i=0; i<N; ++i)
      {
      const uword row = i + row_offset;
      const uword col = i + col_offset;
      
      const eT val = X.at(row,col);
      
      if(val != eT(0))  { out.at(row,col) = val; }
      }
    }
  }



//! @}
