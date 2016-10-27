// Copyright (C) 2009-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by Dimitrios Bouzas



//! \addtogroup op_repmat
//! @{



template<typename obj>
inline
void
op_repmat::apply_noalias(Mat<typename obj::elem_type>& out, const obj& X, const uword copies_per_row, const uword copies_per_col)
  {
  arma_extra_debug_sigprint();
  
  typedef typename obj::elem_type eT;
  
  const uword X_n_rows = obj::is_row ? uword(1) : X.n_rows;
  const uword X_n_cols = obj::is_col ? uword(1) : X.n_cols;
  
  out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
  
  const uword out_n_rows = out.n_rows;
  const uword out_n_cols = out.n_cols;
  
  // if( (out_n_rows > 0) && (out_n_cols > 0) )
  //   {
  //   for(uword col = 0; col < out_n_cols; col += X_n_cols)
  //   for(uword row = 0; row < out_n_rows; row += X_n_rows)
  //     {
  //     out.submat(row, col, row+X_n_rows-1, col+X_n_cols-1) = X;
  //     }
  //   }
  
  if( (out_n_rows > 0) && (out_n_cols > 0) )
    {
    if(copies_per_row != 1)
      {
      for(uword col_copy=0; col_copy < copies_per_col; ++col_copy)
        {
        const uword out_col_offset = X_n_cols * col_copy;
        
        for(uword col=0; col < X_n_cols; ++col)
          {
                eT* out_colptr = out.colptr(col + out_col_offset);
          const eT* X_colptr   = X.colptr(col);
          
          for(uword row_copy=0; row_copy < copies_per_row; ++row_copy)
            {
            const uword out_row_offset = X_n_rows * row_copy;
            
            arrayops::copy( &out_colptr[out_row_offset], X_colptr, X_n_rows );
            }
          }
        }
      }
    else
      {
      for(uword col_copy=0; col_copy < copies_per_col; ++col_copy)
        {
        const uword out_col_offset = X_n_cols * col_copy;
        
        for(uword col=0; col < X_n_cols; ++col)
          {
                eT* out_colptr = out.colptr(col + out_col_offset);
          const eT* X_colptr   = X.colptr(col);
          
          arrayops::copy( out_colptr, X_colptr, X_n_rows );
          }
        }
      }
    }
  
  }



template<typename T1>
inline
void
op_repmat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword copies_per_row = in.aux_uword_a;
  const uword copies_per_col = in.aux_uword_b;
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_repmat::apply_noalias(tmp, U.M, copies_per_row, copies_per_col);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_repmat::apply_noalias(out, U.M, copies_per_row, copies_per_col);
    }
  }



//! @}
