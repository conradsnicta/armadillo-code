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


//! \addtogroup spop_normalise
//! @{



template<typename T1>
inline
void
spop_normalise::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_normalise>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword p   = expr.aux_uword_a;
  const uword dim = expr.aux_uword_b;
  
  arma_debug_check( (p   == 0), "normalise(): parameter 'p' must be greater than zero" );
  arma_debug_check( (dim >  1), "normalise(): parameter 'dim' must be 0 or 1"          );
  
  const unwrap_spmat<T1> U(expr.m);
  
  if(&out == &U.M)
    {
    SpMat<eT> tmp;
    
    spop_normalise::apply_noalias(tmp, U.M, p, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    spop_normalise::apply_noalias(out, U.M, p, dim);
    }
  }



template<typename eT>
inline
void
spop_normalise::apply_noalias(SpMat<eT>& out, const SpMat<eT>& X, const uword p, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: this is a brute force implementation that doesn't explicitly handle potential overflows/underflows;
  // NOTE: however, since we're dealing with a sparse matrix (ie. a relatively small amount of non-zeros),
  // NOTE: overflows/underflows are less likely
  // 
  // NOTE: a seemingly more efficient approach would be via inplace modification of elements via an iterator;
  // NOTE: however, that approach is not robust, as writing a zero through an iterator will invalidate the iterator;
  // NOTE: there are no guarantees that a normalised value will be non-zero, eg. in double precision floats: 1.0 / 1e312 = 0
  
  out.copy_size(X);
  
  typename SpMat<eT>::const_iterator it_begin = X.begin();
  typename SpMat<eT>::const_iterator it_end   = X.end();
  
  typename SpMat<eT>::const_iterator it = it_begin;
  
  if(dim == 0)
    {
    podarray<eT> tmp(out.n_cols);
    tmp.zeros();
    
    eT* tmp_mem = tmp.memptr();
    
    if(p == uword(1))
      {
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.col() ] += std::abs(val); }
      }
    else
    if(p == uword(2))
      {
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.col() ] += (val*val); }
      
      for(uword i=0; i<tmp.n_elem; ++i)  { eT& val = tmp_mem[i]; val = std::sqrt(val); }
      }
    else
      {
      const eT pp = eT(p);
      
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.col() ] += std::pow(val, pp); }
      
      const eT inv_pp = eT(1) / pp;
      
      for(uword i=0; i<tmp.n_elem; ++i)  { eT& val = tmp_mem[i]; val = std::pow(val, inv_pp); }
      }
    
    it = it_begin;
    
    for(; it != it_end; ++it)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      out.at(row, col) = (*it) / tmp_mem[col];
      }
    }
  else
  if(dim == 1)
    {
    podarray<eT> tmp(out.n_rows);
    tmp.zeros();
    
    eT* tmp_mem = tmp.memptr();
    
    if(p == uword(1))
      {
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.row() ] += std::abs(val); }
      }
    else
    if(p == uword(2))
      {
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.row() ] += (val*val); }
      
      for(uword i=0; i<tmp.n_elem; ++i)  { eT& val = tmp_mem[i]; val = std::sqrt(val); }
      }
    else
      {
      const eT pp = eT(p);
      
      for(; it != it_end; ++it)  { const eT val = (*it); tmp_mem[ it.row() ] += std::pow(val, pp); }
      
      const eT inv_pp = eT(1) / pp;
      
      for(uword i=0; i<tmp.n_elem; ++i)  { eT& val = tmp_mem[i]; val = std::pow(val, inv_pp); }
      }
    
    it = it_begin;
    
    for(; it != it_end; ++it)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      out.at(row, col) = (*it) / tmp_mem[row];
      }
    }
  }



//! @}
