// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_cumprod
//! @{



template<typename eT>
inline
void
op_cumprod::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  uword n_rows = X.n_rows;
  uword n_cols = X.n_cols;
  
  out.set_size(n_rows,n_cols);
  
  if(dim == 0)
    {
    if(n_cols == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      eT acc = eT(1);
      
      for(uword row=0; row < n_rows; ++row)
        {
        acc *= X_mem[row];
        
        out_mem[row] = acc;
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        const eT*   X_colmem =   X.colptr(col);
              eT* out_colmem = out.colptr(col);
        
        eT acc = eT(1);
        
        for(uword row=0; row < n_rows; ++row)
          {
          acc *= X_colmem[row];
          
          out_colmem[row] = acc;
          }
        }
      }
    }
  else
  if(dim == 1)
    {
    if(n_rows == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      eT acc = eT(1);
      
      for(uword col=0; col < n_cols; ++col)
        {
        acc *= X_mem[col];
        
        out_mem[col] = acc;
        }
      }
    else
      {
      if(n_cols > 0)
        {
        arrayops::copy( out.colptr(0), X.colptr(0), n_rows );
        
        for(uword col=1; col < n_cols; ++col)
          {
          const eT* out_colmem_prev = out.colptr(col-1);
                eT* out_colmem      = out.colptr(col  );
          const eT*   X_colmem      =   X.colptr(col  );
          
          for(uword row=0; row < n_rows; ++row)
            {
            out_colmem[row] = out_colmem_prev[row] * X_colmem[row];
            }
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_cumprod::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumprod>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_debug_check( (dim > 1), "cumprod(): parameter 'dim' must be 0 or 1" );
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_cumprod::apply_noalias(tmp, U.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_cumprod::apply_noalias(out, U.M, dim);
    }
  }



template<typename T1>
inline
void
op_cumprod_default::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumprod_default>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(in.m);
  
  const uword dim = (T1::is_row) ? 1 : 0;
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_cumprod::apply_noalias(tmp, U.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_cumprod::apply_noalias(out, U.M, dim);
    }
  }



//! @}

